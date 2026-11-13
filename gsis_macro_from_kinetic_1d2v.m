function [macro,aux,hdg,hot] = gsis_macro_from_kinetic_1d2v(kinetic,param)
%GSIS_MACRO_FROM_KINETIC_1D2V One synthetic HDG subsolve from kinetic data.
% kinetic needs mesh.edges, H, velocity, and config.delta,
% config.acceleration. New mono1D2V outputs supply these fields directly.
% Full-state coupling patterned after GSIS_Gray_20261110: weak upwind HOT
% derivatives and affine half-range kinetic wall fluxes, with BGK coefficients.
% An additional seven-moment weak load correction matches the original
% kinetic-DG residual at the frozen H state; this is beyond the Gray port.
% The macro operator/stabilization is retained and H is never changed here.
% param.kineticConsistency=false disables ONLY this additional weak load.
% Does not advance or replace the CIS outer loop in mono1D2V.
if nargin<2; param=struct; end
mesh=kinetic.mesh; v=kinetic.velocity; N=numel(mesh.edges)-1;
assert(size(kinetic.H,1)==2*N,'Inconsistent kinetic spatial DOFs.');
hot=compute_high_order_moments_1d2v(kinetic.H,v);
[trace,wall]=kinetic_face_trace_1d2v(kinetic.H,v,kinetic.config.wallVelocity);
faceHot=compute_high_order_moments_1d2v(trace,v);
hot.derivative=project_hot_derivatives_1d(mesh,hot.R,faceHot.R,'weak');
hot.derivativeMethod='upwindWeak'; hot.faceMoments=faceHot.R;
hot.kineticTrace=trace; hot.wallInflow=wall;
hot.Pi_xx=hot.derivative(:,1); hot.Pi_xy=hot.derivative(:,2);
hot.Theta_x=hot.derivative(:,3);
bc=kinetic_macro_wall_bc_1d2v(kinetic,trace);
if ~isfield(param,'cache'); param.cache=hdg_precompute_p1(mesh,kinetic.config.delta,param); end
enabled=true;
if isfield(param,'kineticConsistency')
    validateattributes(param.kineticConsistency,{'logical','numeric'},{'scalar','binary'});
    enabled=logical(param.kineticConsistency);
end
if enabled
    hot.consistency=kinetic_hdg_consistency_1d2v(kinetic,hot,param.cache,bc);
    hot.consistency.enabled=true;
else
    hot.consistency=struct('enabled',false,'method','disabled', ...
        'loadCorrection',zeros(14,N));
end
param.localLoadCorrection=hot.consistency.loadCorrection;
[macro,aux,hdg]=hdg_solve_macro_1d2v(mesh,kinetic.config.delta, ...
    [0 kinetic.config.acceleration],hot.Pi_xx,hot.Pi_xy,hot.Theta_x,bc,param);
hot.boundary=bc;
hot.velocityMoments=[sum(v.weight),sum(v.weight.*v.vx.^2), ...
    sum(v.weight.*v.vx.^2.*v.vy.^2)];
hot.velocityMomentDefect=hot.velocityMoments-[1 .5 .25];
end
