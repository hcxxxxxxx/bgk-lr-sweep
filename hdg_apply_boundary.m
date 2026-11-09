function [trace,info] = hdg_apply_boundary(K,rhs,macro_bc,context)
%HDG_APPLY_BOUNDARY Apply primitive traces or prescribed outward HDG flux.
% macro_bc may be a struct or @(context) returning that struct.
% .trace is 2x4 [p ux uy T], p=rho+T (left/right rows); NaN means unspecified.
% .flux is 2x4 outward [mass 2*normalMomentum 2*shearMomentum heat] flux.
% Each boundary component must prescribe exactly one of trace and flux.
% Flux-only/gauge-free data may be singular; the caller must close the gauge.
if isa(macro_bc,'function_handle'); macro_bc=macro_bc(context); end
assert(isstruct(macro_bc) && isscalar(macro_bc),'Invalid macro boundary data.');
traceData=nan(2,4); fluxData=nan(2,4);
if isfield(macro_bc,'trace'); traceData=macro_bc.trace; end
if isfield(macro_bc,'flux'); fluxData=macro_bc.flux; end
for value={traceData,fluxData}
    assert(isnumeric(value{1}) && isreal(value{1}) && isequal(size(value{1}),[2 4]) ...
        && all(isfinite(value{1}) | isnan(value{1}),'all'), ...
        'Boundary trace/flux arrays must be real 2x4 arrays (NaN is allowed).');
end
assert(all(xor(isfinite(traceData),isfinite(fluxData)),'all'), ...
    'Prescribe exactly one trace or outward flux per boundary component.');
boundaryIds=[(1:4)';(size(K,1)-3:size(K,1))'];
prescribed=reshape(traceData',[],1); targetFlux=reshape(fluxData',[],1);
known=boundaryIds(isfinite(prescribed)); knownValues=prescribed(isfinite(prescribed));
natural=boundaryIds(isfinite(targetFlux));
rhs(natural)=rhs(natural)+targetFlux(isfinite(targetFlux));
free=setdiff((1:size(K,1))',known);
trace=zeros(size(rhs)); trace(known)=knownValues;
if ~isempty(free)
    trace(free)=K(free,free)\(rhs(free)-K(free,known)*knownValues);
end
assert(all(isfinite(trace)),'HDG trace solve failed; check boundary constraints/gauges.');
info=struct('known',known,'knownValues',knownValues,'free',free, ...
    'rhs',rhs,'trace',traceData,'flux',fluxData, ...
    'residual',max([0;abs(K(free,:)*trace-rhs(free));abs(trace(known)-knownValues)]));
end
