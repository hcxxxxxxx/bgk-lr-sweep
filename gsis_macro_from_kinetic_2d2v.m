function [macro,info]=gsis_macro_from_kinetic_2d2v(kinetic,param)
% Full-state 2D2V BGK Q1-HDG subsolve with kinetic affine diffuse-wall fluxes.
% Frozen HOT contains volume AND kinetic upwind face moments. No H lift.
% With transportWallDensity supplied, HOT uses the actual last-sweep wall
% inflow, while the affine macro wall flux uses newly reflected post-sweep H.
% Default: no added internal residual load. param.kineticConsistency=true
% optionally matches all eight kinetic weak moment residuals at current H.
if nargin<2;param=struct;end
mesh=kinetic.mesh;v=kinetic.velocity;H=kinetic.H;delta=kinetic.config.delta;N=mesh.count;
if isfield(param,'cache');c=param.cache;else;c=hdg_precompute_q1_2d2v(mesh,v,delta,param);end
assert(c.delta==delta && isequal(c.edgesX,mesh.xEdges) && isequal(c.edgesY,mesh.yEdges));
assert(isfield(c,'formulation') && strcmp(c.formulation,'gray-bgk-2v-v2'),'Rebuild the 2D macro cache for the current formulation.');
if isfield(param,'stabilization');assert(isequal(c.stabilization,param.stabilization(:)'));end
q=H*v.moment;primitive=[q(:,1)+q(:,4),q(:,2:4)];
x=v.vx;y=v.vy;s=v.speed2;
hotBasis=[x.*(x.^2-y.^2-1),y.*(x.^2-y.^2+1), ...
    y.*(2*x.^2-1),x.*(2*y.^2-1), ...
    x.^2.*(s-2)-(s-1),x.*y.*(s-2), ...
    x.*y.*(s-2),y.^2.*(s-2)-(s-1)];
fluxBasis=[ones(size(x)),2*x,2*y,s-2];R=H*(v.weight.*hotBasis);
enabled=false;
if isfield(param,'kineticConsistency')
    validateattributes(param.kineticConsistency,{'logical','numeric'},{'scalar','binary'});
    enabled=logical(param.kineticConsistency);
end
test=[fluxBasis,x.^2-y.^2,2*x.*y,x.*(s-2),y.*(s-2)];
momentResidual=zeros(32,N);baseFields=zeros(32,N);wallFlux=cell(N,4);
if enabled
    gain=(q/v.gram)*v.basis';source=delta*(gain-H)*(v.weight.*test);
    fluxX=H*((v.weight.*x).*test);fluxY=H*((v.weight.*y).*test);
    stress=H*(v.weight.*test(:,5:8));
end
localNodes={[1 3],[2 4],[1 2],[3 4]};opposite=[2 1 4 3];
localRhs=zeros(32,N);base=zeros(32,N);rhs=zeros(8*c.faceCount,1);meanBase=0;
wallTargets=cell(N,4);hotLoad=zeros(4,4,N);
hasSweepWall=isfield(kinetic,'transportWallDensity');
if hasSweepWall;assert(isequal(size(kinetic.transportWallDensity),[4*N 4]));end
for j=1:N
    ids=4*(j-1)+(1:4);load=-c.Gx(:,:,j)*R(ids,[1 3 5 7])-c.Gy(:,:,j)*R(ids,[2 4 6 8]);
    if enabled
        kineticResidual=-c.Gx(:,:,j)*fluxX(ids,:)-c.Gy(:,:,j)*fluxY(ids,:)-c.M(:,:,j)*source(ids,:);
        baseFields(:,j)=reshape([primitive(ids,:),-delta*stress(ids,:)],[],1);
    end
    for side=1:4
        normal=c.normals(side,:);cn=normal(1)*x+normal(2)*y;incoming=cn<0;
        nodes=4*(j-1)+localNodes{side};trace=H(nodes,:);neighbor=mesh.neighbor(j,side);
        if neighbor>0
            external=4*(neighbor-1)+localNodes{opposite(side)};
            trace(:,incoming)=H(external,incoming);
        else
            outgoing=cn>0;motion=zeros(1,numel(x));
            if side==4;motion=2*kinetic.config.lidSpeed*x';end
            rho=-(trace(:,outgoing)*(v.weight(outgoing).*cn(outgoing)) ...
                +sum(motion(incoming)'.*v.weight(incoming).*cn(incoming)))/sum(v.weight(incoming).*cn(incoming));
            trace(:,incoming)=rho+motion(incoming);
            F=trace*((v.weight.*cn).*fluxBasis);
            target=c.faceM(:,:,side,j)*(F-primitive(nodes,:)*c.wallJ(:,:,side)');
            wallTargets{j,side}=target;
            wallFlux{j,side}=c.faceM(:,:,side,j)*F;
        end
        hotTrace=trace;
        if neighbor==0 && hasSweepWall
            sweepDensity=kinetic.transportWallDensity(nodes,side);
            assert(all(isfinite(sweepDensity)),'Missing last-sweep incoming wall density.');
            hotTrace(:,incoming)=sweepDensity+motion(incoming);
        end
        faceR=hotTrace*(v.weight.*hotBasis);
        normalR=normal(1)*faceR(:,[1 3 5 7])+normal(2)*faceR(:,[2 4 6 8]);
        load=load+c.E{side}*c.faceM(:,:,side,j)*normalR;
        if enabled
            kineticResidual=kineticResidual+c.E{side}*c.faceM(:,:,side,j)*(trace*((v.weight.*cn).*test));
        end
    end
    localRhs(17:32,j)=load(:);hotLoad(:,:,j)=load;
    if enabled;momentResidual(:,j)=reshape([kineticResidual(:,1:4),-kineticResidual(:,5:8)],[],1);end
    base(:,j)=c.A(:,:,j)\localRhs(:,j);
    density=zeros(1,32);density(1:4)=sum(c.M(:,:,j),1);density(13:16)=-sum(c.M(:,:,j),1);
    meanBase=meanBase+density*base(:,j);
    for side=1:4
        localRows=8*(side-1)+(1:8);globalRows=c.ids(localRows,j);
        rhs(globalRows)=rhs(globalRows)-c.CJ(:,:,side,j)*base(:,j);
        if mesh.neighbor(j,side)==0;rhs(globalRows)=rhs(globalRows)+wallTargets{j,side}(:);end
    end
end
loadCorrection=zeros(32,N);baseTrace=[];
if enabled
    % Choose a base trace satisfying HDG face conservation and the actual
    % kinetic wall flux. Freeze the remaining local discrepancy as a load.
    faceD=zeros(8,8,c.faceCount);faceFlux=zeros(8,c.faceCount);baseTrace=zeros(8,c.faceCount);
    for j=1:N
        for side=1:4
            face=c.faces(j,side);faceD(:,:,face)=faceD(:,:,face)+c.D(:,:,side,j);
            faceFlux(:,face)=faceFlux(:,face)+c.C(:,:,side,j)*baseFields(:,j);
        end
    end
    for face=1:c.faceCount;baseTrace(:,face)=-faceD(:,:,face)\faceFlux(:,face);end
    for j=1:N
        for side=1:4
            if mesh.neighbor(j,side)>0;continue;end
            face=c.faces(j,side);
            baseTrace(:,face)=c.D(:,:,side,j)\(wallFlux{j,side}(:)-c.C(:,:,side,j)*baseFields(:,j));
        end
    end
    rhs(:)=0;meanBase=0;
    for j=1:N
        hat0=reshape(baseTrace(:,c.faces(j,:)),[],1);
        loadCorrection(:,j)=c.A(:,:,j)*baseFields(:,j)-c.B(:,:,j)*hat0-momentResidual(:,j)-localRhs(:,j);
        localRhs(:,j)=localRhs(:,j)+loadCorrection(:,j);
        base(:,j)=c.A(:,:,j)\localRhs(:,j);
        density=zeros(1,32);density(1:4)=sum(c.M(:,:,j),1);density(13:16)=-sum(c.M(:,:,j),1);
        meanBase=meanBase+density*base(:,j);
        for side=1:4
            localRows=8*(side-1)+(1:8);globalRows=c.ids(localRows,j);
            rhs(globalRows)=rhs(globalRows)-c.CJ(:,:,side,j)*base(:,j);
            if mesh.neighbor(j,side)==0;rhs(globalRows)=rhs(globalRows)+wallTargets{j,side}(:);end
        end
    end
end
sol=c.factor\[rhs;-meanBase];hat=sol(1:end-1);Z=zeros(32,N);localResidual=0;
flux=zeros(8*c.faceCount,1);wallResidual=0;
for j=1:N
    h=hat(c.ids(:,j));Z(:,j)=base(:,j)+c.response(:,:,j)*h;
    localResidual=max(localResidual,max(abs(c.A(:,:,j)*Z(:,j)-localRhs(:,j)-c.B(:,:,j)*h)));
    for side=1:4
        localRows=8*(side-1)+(1:8);globalRows=c.ids(localRows,j);
        F=c.CJ(:,:,side,j)*Z(:,j)+c.D(:,:,side,j)*h(localRows);
        if mesh.neighbor(j,side)==0
            F=F-wallTargets{j,side}(:);wallResidual=max(wallResidual,max(abs(F)));
        end
        flux(globalRows)=flux(globalRows)+F;
    end
end
P=reshape(permute(reshape(Z(1:16,:),4,4,N),[1 3 2]),4*N,4);
aux=reshape(permute(reshape(Z(17:32,:),4,4,N),[1 3 2]),4*N,4);
macro=struct('pressureNodal',P,'nodal',[P(:,1)-P(:,4),P(:,2:4)],'auxiliary',aux);
macro.cellAverage=reshape(mean(reshape(macro.nodal,4,N,4),1),N,4);
info=struct('localResidual',localResidual,'traceResidual',max(abs(flux)), ...
    'wallFluxResidual',wallResidual,'densityGaugeMultiplier',sol(end), ...
    'linearResidual',max(abs(c.system*sol-[rhs;-meanBase])), ...
    'localFields',Z,'trace',hat,'localRhs',localRhs,'hotWeakLoad',hotLoad, ...
    'formulation',c.formulation,'stabilization',c.stabilization,'internalResidualCorrection',enabled, ...
    'hotUsesTransportWallInflow',hasSweepWall, ...
    'loadCorrection',loadCorrection,'baseFields',baseFields,'baseTrace',baseTrace, ...
    'targetKineticMomentResidual',momentResidual);
info.wallAffineTarget=wallTargets;info.wallKineticFlux=wallFlux;
end
