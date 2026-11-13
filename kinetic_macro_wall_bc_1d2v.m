function bc = kinetic_macro_wall_bc_1d2v(kinetic,trace)
%KINETIC_MACRO_WALL_BC_1D2V F_HDG=F_DVM+J*(U-U_DVM) for diffuse walls.
% U=[p ux uy T], p=rho+T. Flux rows are [mass,2*normalMomentum,
% 2*shearMomentum,heat]. The tangent uses the SAME discrete Gram projection
% as BGK and the SAME discrete zero-normal-mass diffuse law as kinetic.
v=kinetic.velocity; H=kinetic.H;
if nargin<2
    trace=kinetic_face_trace_1d2v(H,v,kinetic.config.wallVelocity);
end
q=H([1 end],:)*v.moment;
state=[q(:,1)+q(:,4),q(:,2:4)];
fluxBasis=[ones(size(v.vx)),2*v.vx,2*v.vy,v.speed2-2];
% Row U -> row [rho ux uy T] -> velocity gain.
primitiveToMoments=[1 0 0 0;0 1 0 0;0 0 1 0;-1 0 0 1];
tangent=(primitiveToMoments/v.gram*v.basis')';
flux=zeros(2,4); jacobian=zeros(4,4,2);
for side=1:2
    cn=(2*side-3)*v.vx; outgoing=cn>0; incoming=cn<0;
    wallTangent=tangent;
    reflected=-sum((v.weight(outgoing).*cn(outgoing)).*tangent(outgoing,:),1) ...
        /sum(v.weight(incoming).*cn(incoming));
    wallTangent(incoming,:)=repmat(reflected,nnz(incoming),1);
    weightedFlux=(v.weight.*cn).*fluxBasis;
    jacobian(:,:,side)=weightedFlux'*wallTangent;
    flux(side,:)=trace(1+(side-1)*(size(trace,1)-1),:)*weightedFlux;
end
bc=struct('kind','kineticFlux','flux',flux,'state',state,'jacobian',jacobian, ...
    'meanDensity',0);
end
