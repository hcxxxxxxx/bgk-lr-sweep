function bc = macro_bc_from_kinetic(boundaryTrace,velocity)
%MACRO_BC_FROM_KINETIC Reference-code wall data from the INTERIOR distribution.
% boundaryTrace: H([1 end],:), not the upwind/diffuse numerical trace.
% Pressure/temperature and tangential slip are retained; wall-normal ux=0.
% Returns HDG boundary columns [p ux uy T], with p=rho+T.
vx=velocity.vx(:); vy=velocity.vy(:); w=velocity.weight(:);
assert(isequal(size(boundaryTrace),[2 numel(w)]) ...
    && isreal(boundaryTrace) && all(isfinite(boundaryTrace),'all'), ...
    'Provide left/right interior boundary distributions, 2-by-Nv.');
basis=[ones(numel(w),1),vx,vy,vx.^2+vy.^2-1];
moments=boundaryTrace*(w.*basis);
bc=struct('trace',[moments(:,1)+moments(:,4),zeros(2,1),moments(:,3:4)]);
end
