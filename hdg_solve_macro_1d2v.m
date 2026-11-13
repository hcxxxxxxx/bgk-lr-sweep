function [macro,aux,hdg] = hdg_solve_macro_1d2v(mesh,delta,acceleration,Pi_xx,Pi_xy,Theta_x,macro_bc,param)
%HDG_SOLVE_MACRO_1D2V Synthetic 1D2V BGK moment equations, P1 HDG.
% [macro,aux,hdg]=hdg_solve_macro_1d2v(mesh,delta,[ax ay], ...
%     Pi_xx,Pi_xy,Theta_x,macro_bc,param)
% Corrections: scalar, 2*N nodal column, N-by-2 cell/node array, or @(x).
% Acceleration: [ax ay], a 2*N-by-2 nodal array, or @(x) returning that array.
% Boundary interface: hdg_apply_boundary for trace/prescribed flux data;
% kind='kineticFlux' uses affine flux/Jacobian data and a mean-density gauge.
% Boundary trace columns and local fields use [p ux uy T], p=rho+T.
% param.stabilization=[s_p s_ux s_uy s_T] defaults to [1 delta delta delta].
% param.cache can reuse hdg.cache for unchanged edges, delta and stabilization.
% param.localLoadCorrection optionally adds a frozen 14-by-N weak load.
% Outputs use 2*N nodal columns, matching mono1D2V; Uhat is (N+1)-by-4.
% sigma_xy = 2*Pxy; aux.Pxy=-Ly/(2*delta), aux.qx=-E/delta.
% This is the macroscopic subsolve, not a kinetic/GSIS outer iteration.
if nargin<8; param=struct; end
if isfield(param,'cache')
    cache=param.cache; s=cache.stabilization;
    if isfield(param,'stabilization'); s=param.stabilization(:)'; end
    assert(isequal(cache.edges,mesh.edges(:)) && isequal(cache.delta,delta) ...
        && isequal(cache.stabilization,s) && isfield(cache,'formulation') ...
        && strcmp(cache.formulation,'pressure-reference-v1'), ...
        'HDG cache does not match the problem/formulation.');
else
    cache=hdg_precompute_p1(mesh,delta,param);
end
N=cache.count; x=cache.nodes;
loadCorrection=zeros(14,N);
if isfield(param,'localLoadCorrection')
    loadCorrection=param.localLoadCorrection;
    assert(isequal(size(loadCorrection),[14 N]) && all(isfinite(loadCorrection),'all'), ...
        'localLoadCorrection must be a finite 14-by-N weak load.');
end
correction=[nodal_field(Pi_xx,N,x),nodal_field(Pi_xy,N,x),nodal_field(Theta_x,N,x)];
if isa(acceleration,'function_handle'); acceleration=acceleration(x); end
if isvector(acceleration) && numel(acceleration)==2
    acceleration=repmat(acceleration(:)',2*N,1);
end
assert(isnumeric(acceleration) && isreal(acceleration) ...
    && isequal(size(acceleration),[2*N 2]) && all(isfinite(acceleration),'all'), ...
    'acceleration must specify [ax ay] or a 2*N-by-2 nodal array.');
forcing=[zeros(2*N,1),2*acceleration,zeros(2*N,1),correction];
localRhs=zeros(14,N); U0=zeros(14,N); rhs=zeros(4*(N+1),1);
for i=1:N
    nodes=2*i-1:2*i;
    localRhs(:,i)=reshape(cache.M(:,:,i)*forcing(nodes,:),[],1)+loadCorrection(:,i);
    U0(:,i)=cache.A(:,:,i)\localRhs(:,i);
    ids=cache.traceIds(:,i);
    rhs(ids)=rhs(ids)-cache.C(:,:,i)*U0(:,i);
end
context=struct('mesh',mesh,'delta',delta,'acceleration',acceleration, ...
    'correction',correction,'boundaryX',cache.edges([1 end]),'normals',[-1;1]);
if isa(macro_bc,'function_handle'); macro_bc=macro_bc(context); end
kineticBoundary=isfield(macro_bc,'kind') && strcmp(macro_bc.kind,'kineticFlux');
if kineticBoundary
    [trace,bc]=hdg_apply_kinetic_boundary_1d2v(cache,U0,rhs,macro_bc);
else
    [trace,bc]=hdg_apply_boundary(cache.K,rhs,macro_bc,context);
end
U=zeros(14,N); outwardFlux=zeros(8,N); localResidual=0;
for i=1:N
    hat=trace(cache.traceIds(:,i));
    U(:,i)=U0(:,i)+cache.response(:,:,i)*hat;
    outwardFlux(:,i)=cache.C(:,:,i)*U(:,i)+cache.D(:,:,i)*hat;
    localResidual=max(localResidual,max(abs(cache.A(:,:,i)*U(:,i)-localRhs(:,i)-cache.B(:,:,i)*hat)));
end
macro=struct; names={'p','ux','uy','T'};
for k=1:4; macro.(names{k})=reshape(U(2*k-1:2*k,:),[],1); end
macro.rho=macro.p-macro.T;
macro.pressureNodal=[macro.p,macro.ux,macro.uy,macro.T];
macro.nodal=[macro.rho,macro.ux,macro.uy,macro.T];
macro.cellAverage=(macro.nodal(1:2:end,:)+macro.nodal(2:2:end,:))/2;
aux=struct('Lx',reshape(U(9:10,:),[],1),'Ly',reshape(U(11:12,:),[],1), ...
    'E',reshape(U(13:14,:),[],1));
aux.sigma_xx=-aux.Lx/delta; aux.sigma_xy=-aux.Ly/delta;
aux.Pxy=aux.sigma_xy/2; aux.qx=-aux.E/delta;
fluxJump=outwardFlux(5:8,1:end-1)+outwardFlux(1:4,2:end);
cellBalance=zeros(N,4);
for i=1:N
    source=cache.M(:,:,i)*forcing(2*i-1:2*i,1:4)+reshape(loadCorrection(1:8,i),2,4);
    cellBalance(i,:)=(outwardFlux(1:4,i)+outwardFlux(5:8,i))'-sum(source,1);
end
hdg=struct('Uhat',reshape(trace,4,[])','localFields',U,'localRhs',localRhs, ...
    'outwardFlux',outwardFlux,'fluxJump',fluxJump,'cellBalance',cellBalance, ...
    'traceResidual',bc.residual,'localResidual',localResidual, ...
    'fluxContinuityResidual',max([0;abs(fluxJump(:))]), ...
    'conservationResidual',max(abs(cellBalance),[],'all'), ...
    'cache',cache,'bc',bc,'traceRhs',rhs,'U0',U0);
hdg.traceColumns={'p','ux','uy','T'};
if kineticBoundary
    mismatch=zeros(2,4);
    for side=1:2
        cellId=1+(side-1)*(N-1); localRows=(side-1)*4+(1:4);
        primitive=U(side:2:8,cellId);
        expected=macro_bc.flux(side,:)'+macro_bc.jacobian(:,:,side)*(primitive-macro_bc.state(side,:)');
        mismatch(side,:)=outwardFlux(localRows,cellId)'-expected';
    end
    hdg.boundaryFluxResidual=max(abs(mismatch),[],'all');
    hdg.boundaryFluxMismatch=mismatch;
end
end

function value=nodal_field(value,N,x)
if isa(value,'function_handle'); value=value(x); end
if isscalar(value)
    value=repmat(value,2*N,1);
elseif isequal(size(value),[N 2])
    value=reshape(value',[],1);
else
    value=value(:);
end
assert(isnumeric(value) && isreal(value) && numel(value)==2*N ...
    && all(isfinite(value)),'Correction fields must have 2*N nodal entries.');
end
