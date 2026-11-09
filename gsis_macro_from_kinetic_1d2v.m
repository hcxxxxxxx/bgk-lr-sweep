function [macro,aux,hdg,hot] = gsis_macro_from_kinetic_1d2v(kinetic,param)
%GSIS_MACRO_FROM_KINETIC_1D2V One synthetic HDG subsolve from kinetic data.
% kinetic needs mesh.edges, H, velocity, and config.delta,
% config.acceleration. New mono1D2V outputs supply these fields directly.
% Uses element-interior HOT derivatives and interior wall moments, matching
% NonOscillatory/Synthetic_Acceleration.f90, with 1D2V BGK coefficients.
% Does not advance or replace the CIS outer loop in mono1D2V.
if nargin<2; param=struct; end
mesh=kinetic.mesh; v=kinetic.velocity; N=numel(mesh.edges)-1;
assert(size(kinetic.H,1)==2*N,'Inconsistent kinetic spatial DOFs.');
hot=compute_high_order_moments_1d2v(kinetic.H,v);
hot.derivative=project_hot_derivatives_1d(mesh,hot.R,[],'strong');
hot.derivativeMethod='interior';
hot.Pi_xx=hot.derivative(:,1); hot.Pi_xy=hot.derivative(:,2);
hot.Theta_x=hot.derivative(:,3);
boundaryTrace=kinetic.H([1 end],:);
bc=@(~) macro_bc_from_kinetic(boundaryTrace,v);
[macro,aux,hdg]=hdg_solve_macro_1d2v(mesh,kinetic.config.delta, ...
    [0 kinetic.config.acceleration],hot.Pi_xx,hot.Pi_xy,hot.Theta_x,bc,param);
hot.boundary=macro_bc_from_kinetic(boundaryTrace,v);
hot.velocityMoments=[sum(v.weight),sum(v.weight.*v.vx.^2), ...
    sum(v.weight.*v.vx.^2.*v.vy.^2)];
hot.velocityMomentDefect=hot.velocityMoments-[1 .5 .25];
end
