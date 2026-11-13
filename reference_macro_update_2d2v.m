function [residual,component]=reference_macro_update_2d2v(Q,previous,mesh,tol)
% Gray Update_Macro_Residual, adapted to [rho,ux,uy,T] in 2V BGK.
% Use cell INTEGRALS of full density/temperature and velocity; the velocity
% reference scale is one under the present BGK thermal-speed normalization.
% As in the archive, compare each component separately and use the same
% tol-dependent denominator floor. Q is the returned source for GSIS.
N=mesh.count;
current=reshape(mean(reshape(Q,4,N,4),1),N,4);
old=reshape(mean(reshape(previous,4,N,4),1),N,4);
current(:,[1 4])=current(:,[1 4])+1;
old(:,[1 4])=old(:,[1 4])+1;
current=mesh.area.*current;old=mesh.area.*old;
numerator=sum((current-old).^2,1);denominator=sum(current.^2,1);
reference=max([denominator,1]);active=denominator>tol^2*reference;
component=sqrt(numerator/reference);
component(active)=sqrt(numerator(active)./denominator(active));
residual=max(component);
end
