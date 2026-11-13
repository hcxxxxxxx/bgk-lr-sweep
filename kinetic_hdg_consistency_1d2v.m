function consistency = kinetic_hdg_consistency_1d2v(kinetic,hot,cache,bc)
%KINETIC_HDG_CONSISTENCY_1D2V Match seven HDG weak residuals to kinetic moments.
% Freeze discrepancies between stabilized HDG faces and kinetic upwind faces
% in the macro load; retain the HDG operator and its stabilization parameters.
% At the kinetic state the conservation residuals equal four kinetic moments;
% the auxiliary residuals equal minus three kinetic stress/heat moments.
% A converged CIS state therefore receives zero macroscopic correction.
H=kinetic.H; v=kinetic.velocity; N=cache.count; delta=cache.delta;
q=H*v.moment; pressure=[q(:,1)+q(:,4),q(:,2:4)];
test=[ones(size(v.vx)),2*v.vx,2*v.vy,v.speed2-2, ...
    v.vx.^2-v.vy.^2,2*v.vx.*v.vy,v.vx.*(v.speed2-2)];
weighted=v.weight.*test; stress=H*weighted(:,5:7);
baseFields=zeros(14,N);
for cellId=1:N
    nodes=2*cellId-1:2*cellId;
    baseFields(:,cellId)=reshape([pressure(nodes,:),-delta*stress(nodes,:)],[],1);
end
baseTrace=zeros(N+1,4);
for face=2:N
    left=cache.C(5:8,:,face-1)*baseFields(:,face-1);
    right=cache.C(1:4,:,face)*baseFields(:,face);
    baseTrace(face,:)=((left+right)./(2*cache.stabilization(:)))';
end
for side=1:2
    cellId=1+(side-1)*(N-1); rows=(side-1)*4+(1:4);
    baseTrace(1+(side-1)*N,:)=(cache.D(rows,rows,cellId)\ ...
        (bc.flux(side,:)'-cache.C(rows,:,cellId)*baseFields(:,cellId)))';
end
gain=(q/v.gram)*v.basis'; force=2*kinetic.config.acceleration*v.vy';
source=(delta*(gain-H)+force)*weighted;
volumeFlux=H*((v.weight.*v.vx).*test);
faceFlux=hot.kineticTrace*((v.weight.*v.vx).*test);
forcing=[zeros(2*N,2),repmat(2*kinetic.config.acceleration,2*N,1), ...
    zeros(2*N,1),hot.derivative];
loadCorrection=zeros(14,N); kineticResidual=zeros(2*N,7);
for cellId=1:N
    nodes=2*cellId-1:2*cellId; M=cache.M(:,:,cellId);
    residual=-cache.G*volumeFlux(nodes,:)+[-faceFlux(cellId,:);faceFlux(cellId+1,:)] ...
        -M*source(nodes,:);
    kineticResidual(nodes,:)=residual;
    target=reshape([residual(:,1:4),-residual(:,5:7)],[],1);
    hat=reshape(baseTrace(cellId:cellId+1,:)',[],1);
    originalLoad=reshape(M*forcing(nodes,:),[],1);
    loadCorrection(:,cellId)=cache.A(:,:,cellId)*baseFields(:,cellId) ...
        -cache.B(:,:,cellId)*hat-originalLoad-target;
end
consistency=struct('loadCorrection',loadCorrection,'baseFields',baseFields, ...
    'baseTrace',baseTrace,'kineticWeakMomentResidual',kineticResidual, ...
    'method','sevenMomentWeakResidualMatch');
end
