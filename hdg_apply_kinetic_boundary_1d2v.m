function [trace,info] = hdg_apply_kinetic_boundary_1d2v(cache,cellPart,rhs,bc)
%HDG_APPLY_KINETIC_BOUNDARY_1D2V Condense affine kinetic-flux wall conditions.
% All conservation rows are retained. A bordered volume-mean density
% constraint fixes the closed-channel pressure/density null mode.
N=cache.count; K=cache.K;
assert(isequal(size(bc.flux),[2 4]) && isequal(size(bc.state),[2 4]) ...
    && isequal(size(bc.jacobian),[4 4 2]));
for side=1:2
    cellId=1+(side-1)*(N-1); localRows=(side-1)*4+(1:4);
    faceIds=cache.traceIds(localRows,cellId);
    select=zeros(4,14);
    for component=1:4; select(component,2*(component-1)+side)=1; end
    C=cache.C(localRows,:,cellId)-bc.jacobian(:,:,side)*select;
    K(faceIds,:)=0;
    K(faceIds,cache.traceIds(:,cellId))=C*cache.response(:,:,cellId)+cache.D(localRows,:,cellId);
    rhs(faceIds)=bc.flux(side,:)'-bc.jacobian(:,:,side)*bc.state(side,:)' ...
        -C*cellPart(:,cellId);
end
meanRow=zeros(1,size(K,1)); meanBase=0; gap=cache.edges(end)-cache.edges(1);
for cellId=1:N
    measure=sum(cache.M(:,:,cellId),1)/gap;
    density=zeros(1,14); density(1:2)=measure; density(7:8)=-measure;
    ids=cache.traceIds(:,cellId);
    meanRow(ids)=meanRow(ids)+density*cache.response(:,:,cellId);
    meanBase=meanBase+density*cellPart(:,cellId);
end
massCompatibility=repmat([1;0;0;0],N+1,1);
massCompatibility=massCompatibility/norm(massCompatibility);
system=[K,sparse(massCompatibility);sparse(meanRow),sparse(1,1)];
solution=system\[rhs;bc.meanDensity-meanBase]; trace=solution(1:end-1);
assert(all(isfinite(solution)),'Kinetic-flux HDG solve failed.');
info=struct('kind','kineticFlux','flux',bc.flux,'state',bc.state, ...
    'jacobian',bc.jacobian,'meanDensity',bc.meanDensity, ...
    'gaugeMultiplier',solution(end),'densityConstraintResidual', ...
    abs(meanRow*trace+meanBase-bc.meanDensity),'rhs',rhs, ...
    'residual',max(abs(K*trace-rhs)));
end
