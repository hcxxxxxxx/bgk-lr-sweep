function r=monoCavity2D2V(options)
%MONOCAVITY2D2V Linearized BGK lid-driven unit square; genuine 2D Q1 DG.
% r=monoCavity2D2V(struct('delta',10,'cells',[24 24]));
% x,y in [0,1], top wall ux=lidSpeed, all other wall speeds zero.
% Four isothermal fully diffuse walls. No acceleration or VDF lift.
% iterationMethod='GSIS' returns full 2D HDG moments to the collision source.
% Velocity weights include exp(-vx^2-vy^2)/pi. H and rho,T are perturbations.
% Spatial operators follow monoBGKDG, extended to both interior directions.
if nargin<1;options=struct;end
cfg=struct('cells',[24 24],'velocityOrder',[32 32],'delta',1,'lidSpeed',.01, ...
    'quadrature','referenceCubic','cutoff',4,'tol',1e-7,'maxIter',10000, ...
    'meshStretch',0,'meshMapping','tanh','minIter',0,'tolCoupled',Inf, ...
    'iterationMethod','CIS','krylovRestart',40, ...
    'macroParam',struct, ...
    'reportEvery',50,'plot',true,'save',true, ...
    'outputDir',fullfile(fileparts(mfilename('fullpath')),'results','cavity_2d2v'));
for f=fieldnames(options)'
    assert(isfield(cfg,f{1}),'Unknown option %s',f{1});cfg.(f{1})=options.(f{1});
end
if isscalar(cfg.cells);cfg.cells=[cfg.cells cfg.cells];end
validateattributes(cfg.cells,{'numeric'},{'vector','numel',2,'integer','>=',2});
validateattributes(cfg.velocityOrder,{'numeric'},{'vector','numel',2,'integer','>=',4});
validateattributes(cfg.delta,{'numeric'},{'scalar','positive','finite'});
validateattributes(cfg.lidSpeed,{'numeric'},{'scalar','real','finite'});
validateattributes(cfg.tol,{'numeric'},{'scalar','positive','finite'});
validateattributes(cfg.maxIter,{'numeric'},{'scalar','integer','positive'});
validateattributes(cfg.minIter,{'numeric'},{'scalar','integer','nonnegative','<=',cfg.maxIter});
validateattributes(cfg.tolCoupled,{'numeric'},{'scalar','real','positive','nonnan'});
cfg.meshMapping=char(validatestring(cfg.meshMapping,{'tanh','quintic','septic','nonic'}));
cfg.quadrature=char(validatestring(cfg.quadrature,{'referenceCubic','hermite'}));
cfg.iterationMethod=char(validatestring(cfg.iterationMethod,{'CIS','GMRES','GSIS'}));
assert(isstruct(cfg.macroParam) && isscalar(cfg.macroParam),'macroParam must be a scalar struct.');
validateattributes(cfg.meshStretch,{'numeric'},{'scalar','nonnegative','finite'});
v=velocities(cfg);mesh=grid(cfg.cells,cfg.meshStretch,cfg.meshMapping);op=operators(mesh,v,cfg.delta);
N=mesh.count;nv=numel(v.weight);H=zeros(4*N,nv);Q=H*v.moment;
useGSIS=strcmp(cfg.iterationMethod,'GSIS');sourceQ=Q;synthetic=struct;
macroParam=cfg.macroParam;
relaxation=1;requireOriginal=false;
if useGSIS
    if isfield(macroParam,'relaxation');relaxation=macroParam.relaxation;end
    validateattributes(relaxation,{'numeric'},{'scalar','positive','<=',1,'finite'});
    requireOriginal=isfield(macroParam,'kineticConsistency') && logical(macroParam.kineticConsistency);
    if isfield(macroParam,'requireOriginalResidual');requireOriginal=logical(macroParam.requireOriginalResidual);end
    macroParam.cache=hdg_precompute_q1_2d2v(mesh,v,cfg.delta,macroParam);
end
history=zeros(cfg.maxIter,6+7*useGSIS);converged=false;firstResidual=[];timer=tic;
fprintf('%s BGK cavity: %dx%d Q1 cells, %dx%d velocities, delta=%g, lid ux=%g\n', ...
    cfg.iterationMethod,cfg.cells,cfg.velocityOrder,cfg.delta,cfg.lidSpeed);
fprintf('step    relative macro update    relative equation residual    wall mass flux\n');
linearInfo=struct('method',cfg.iterationMethod);transportCalls=0;coupledAbsolute=NaN;coupledRelative=NaN;sourceUpdate=NaN;
referenceComponents=nan(1,4);
if strcmp(cfg.iterationMethod,'GMRES')
    % Unknowns are cell moments AND diffuse-wall densities. Eliminate H by
    % the exact transport sweep; solve the original fixed-point equation.
    % No synthetic macro PDE, no modification of the kinetic equation.
    stateSize=16*N+sum(cellfun(@numel,op.nodes));
    [base,Hfirst]=map(zeros(stateSize,1));
    firstWall=walls(Hfirst,mesh,v,op,cfg.lidSpeed);
    firstGain=((Hfirst*v.moment)/v.gram)*v.basis';
    firstR=-cfg.delta*(op.M*firstGain)-firstWall.load;
    for j=1:nv;firstR(:,j)=firstR(:,j)+op.transport{j}*Hfirst(:,j);end
    firstResidual=max(residual_norm(firstR,op,v),eps*max(abs(cfg.lidSpeed),eps));
    restart=cfg.krylovRestart;cycles=ceil(cfg.maxIter/restart);
    [state,flag,algebraicRelative,krylovIter,krylovHistory]=gmres( ...
        @(z)z-(map(z)-base),base,restart,cfg.tol*.01,cycles,[],[],base);
    [newState,H,used]=map(state);Q=H*v.moment;shift=0;
    update=norm(newState-state)/max(norm(newState),eps);
    iter=transportCalls;
    [absolute,relative,wallFlux,wallChange]=final_diagnostics(H,Q,used,0);
    converged=flag==0 && relative<cfg.tol && update<cfg.tol;
    history=[iter,update,absolute,relative,wallFlux,toc(timer)];
    linearInfo=struct('method','GMRES','flag',flag,'algebraicRelative',algebraicRelative, ...
        'iterations',krylovIter,'residualHistory',krylovHistory, ...
        'transportCalls',transportCalls,'restart',restart);
    fprintf('GMRES flag=%d, transport calls=%d, original relative kinetic residual %.5e\n',flag,iter,relative);
else
for iter=1:cfg.maxIter
    previous=Q;oldSource=sourceQ;used=walls(H,mesh,v,op,cfg.lidSpeed);
    gain=(sourceQ/v.gram)*v.basis';rhs=cfg.delta*(op.M*gain)+used.load;
    for j=1:nv;H(:,j)=op.factor{j}\rhs(:,j);end
    assert(all(isfinite(H),'all'),'Cavity:NonFinite','Non-finite kinetic state; iteration stopped.');
    shift=-(op.integral'*(H*v.weight))/sum(v.weight);H=H+shift;
    Q=H*v.moment;
    % Full density/temperature use the unit reference; velocity uses lid scale.
    scale=[max(norm(1+Q(:,1)),1),max(norm(Q(:,2:3),'fro'),abs(cfg.lidSpeed)*1e-12), ...
        max(norm(1+Q(:,4)),1)];
    difference=Q-previous;
    update=max([norm(difference(:,1))/scale(1), ...
        norm(difference(:,2:3),'fro')/max(scale(2),realmin),norm(difference(:,4))/scale(3)]);
    if useGSIS
        kinetic=struct('H',H,'velocity',v,'mesh',mesh,'config',cfg);
        % Gray HOT reuses the incoming wall load frozen before this sweep;
        % its affine macro wall law separately recomputes the current flux.
        kinetic.transportWallDensity=sweep_wall_density(used,v,op,cfg.lidSpeed,shift);
        [macro,hdg]=gsis_macro_from_kinetic_2d2v(kinetic,macroParam);
        sourceQ=Q+relaxation*(macro.nodal-Q);sourceQ(:,1)=sourceQ(:,1)-op.integral'*sourceQ(:,1);
        assert(all(isfinite(sourceQ),'all'),'Cavity:NonFinite','Non-finite macro source; iteration stopped.');
        sourceUpdate=norm(sourceQ-oldSource,'fro')/max(norm(sourceQ,'fro'),eps);
        % Gray stopping uses the recovered/returned macro cell integrals,
        % not post-sweep kinetic moments or a nodal source update norm.
        [update,referenceComponents]=reference_macro_update_2d2v(sourceQ,oldSource,mesh,cfg.tol);
        synthetic=struct('macro',macro,'hdg',hdg,'sourceMoments',sourceQ, ...
            'relaxation',relaxation,'coupling','macroSourceNoLift', ...
            'referenceMacroUpdateComponents',referenceComponents);
    else
        sourceQ=Q;
    end
    current=walls(H,mesh,v,op,cfg.lidSpeed);
    gain=(Q/v.gram)*v.basis';R=-cfg.delta*(op.M*gain)-current.load;
    for j=1:nv;R(:,j)=R(:,j)+op.transport{j}*H(:,j);end
    absolute=residual_norm(R,op,v);
    if isempty(firstResidual);firstResidual=max(absolute,eps*max(abs(cfg.lidSpeed),eps));end
    relative=absolute/firstResidual;
    coupledAbsolute=absolute;coupledRelative=relative;
    if useGSIS
        sourceGain=(sourceQ/v.gram)*v.basis';
        coupledAbsolute=residual_norm(R-cfg.delta*(op.M*(sourceGain-gain)),op,v);
        coupledRelative=coupledAbsolute/firstResidual;
    end
    % Last sweep's incoming distribution, including the same density gauge.
    wallFlux=0;wallChange=0;
    for s=1:4
        incoming=op.cn{s}<0;outgoing=op.cn{s}>0;nodes=op.nodes{s};
        u=used.state{s}+shift;
        flux=H(nodes,outgoing)*(v.weight(outgoing).*op.cn{s}(outgoing)) ...
            +u(:,incoming)*(v.weight(incoming).*op.cn{s}(incoming));
        wallFlux=max(wallFlux,max(abs(flux)));
        wallChange=max(wallChange,max(abs(u(:,incoming)-current.state{s}(:,incoming)),[],'all'));
    end
    history(iter,1:6)=[iter,update,absolute,relative,wallFlux,toc(timer)];
    if useGSIS;history(iter,7:13)=[sourceUpdate,coupledAbsolute,coupledRelative,referenceComponents];end
    if useGSIS && ~requireOriginal
        % Gray macro stopping by default (tolCoupled=Inf). An explicitly
        % finite tolCoupled also requires the returned-source equation.
        converged=update<cfg.tol && coupledRelative<cfg.tolCoupled;
    else
        converged=update<cfg.tol && coupledRelative<cfg.tol ...
            && wallFlux<cfg.tol*max(abs(cfg.lidSpeed),1e-8);
        if useGSIS;converged=converged && relative<cfg.tol;end
    end
    if iter==1 || mod(iter,cfg.reportEvery)==0 || converged || iter==cfg.maxIter
        fprintf('%5d      %.5e                 %.5e                 %.5e\n',iter,update,relative,wallFlux);drawnow;
        if useGSIS;fprintf('           coupled relative %.5e; macro local/wall %.3e/%.3e\n',coupledRelative,hdg.localResidual,hdg.wallFluxResidual);end
    end
    if converged && iter>=cfg.minIter;break;end
end
end
cellMoments=reshape(mean(reshape(Q,4,N,4),1),N,4);
stress=H*(v.weight.*[v.vx.^2-v.vy.^2,v.vx.*v.vy]);
r=struct('config',cfg,'mesh',mesh,'velocity',v,'H',H,'nodalMoments',Q, ...
    'cellMoments',cellMoments,'shearStress',reshape(mean(reshape(stress(:,2),4,N),1),N,1), ...
    'history',history(1:min(iter,size(history,1)),:),'historyColumns',{{'iteration','macroRelativeUpdate', ...
    'absoluteKineticResidual','relativeKineticResidual','lastSweepWallMassFlux','seconds'}}, ...
    'iterations',iter,'converged',converged,'wallTime',toc(timer));
r.transportWallDensity=sweep_wall_density(used,v,op,cfg.lidSpeed,shift);
r.linearSolver=linearInfo;
if useGSIS
    r.synthetic=synthetic;
    r.historyColumns=[r.historyColumns,{'sourceRelativeUpdate','absoluteCoupledResidual','relativeCoupledResidual', ...
        'referenceDensityUpdate','referenceUxUpdate','referenceUyUpdate','referenceTemperatureUpdate'}];
end
r.diagnostics=struct('kinetic',absolute,'relativeKinetic',relative,'firstResidual',firstResidual, ...
    'macroUpdate',update,'wallMassFlux',wallFlux,'wallInflowChange',wallChange, ...
    'meanDensityPerturbation',op.integral'*Q(:,1), ...
    'collisionInvariantDefect',max(abs((H-gain)*v.moment),[],'all'), ...
    'currentWallMassFlux',current.massFlux);
r.diagnostics.wallMomentumFlux=current.momentumFlux;
r.diagnostics.globalMomentumBalance=max(abs(sum(current.momentumFlux,1)));
r.diagnostics.coupled=coupledAbsolute;r.diagnostics.relativeCoupled=coupledRelative;
r.diagnostics.originalConverged=relative<cfg.tol;
if useGSIS;r.convergenceCriterion='referenceCellIntegratedMacroUpdate';else;r.convergenceCriterion='macroUpdateAndOriginalEquation';end
if useGSIS && requireOriginal;r.convergenceCriterion='referenceMacroUpdateAndOriginalAndCoupledEquation';end
if useGSIS && ~requireOriginal && isfinite(cfg.tolCoupled)
    r.convergenceCriterion='referenceMacroUpdateAndCoupledEquation';
end
% These arrays permit independent small-grid monolithic residual checks.
if cfg.save
    suffix='';
    if ismember(cfg.meshMapping,{'quintic','septic','nonic'});suffix=['_',cfg.meshMapping];
    elseif cfg.meshStretch>0;suffix=sprintf('_stretch%g',cfg.meshStretch);end
    if useGSIS;suffix=[suffix,'_GSIS'];end
    folder=fullfile(cfg.outputDir,[sprintf('delta_%g_N%dx%d',cfg.delta,cfg.cells),suffix]);
    if ~exist(folder,'dir');mkdir(folder);end
    save(fullfile(folder,'cavity.mat'),'r','-v7.3');
    writetable(array2table([mesh.center,cellMoments,r.shearStress], ...
        'VariableNames',{'x','y','rhoPert','ux','uy','TPert','Pxy'}),fullfile(folder,'profile.csv'));
    writetable(array2table(r.history,'VariableNames',r.historyColumns),fullfile(folder,'history.csv'));
    if cfg.plot;plotCavity2D2V(r,folder);end
elseif cfg.plot
    plotCavity2D2V(r);
end
fprintf('Cavity converged=%d, steps=%d, %.2fs, absolute residual %.4e\n', ...
    converged,iter,r.wallTime,absolute);
if useGSIS
    if requireOriginal
        fprintf('GSIS stopping checks both original and coupled equations; original residual criterion passed=%d.\n',r.diagnostics.originalConverged);
    elseif isfinite(cfg.tolCoupled)
        fprintf('GSIS stopping checks macro updates and coupled residual < %.3e; original criterion passed=%d.\n', ...
            cfg.tolCoupled,r.diagnostics.originalConverged);
    else
        fprintf('GSIS stopping uses reference cell-integrated macro updates only; original kinetic residual criterion passed=%d.\n',r.diagnostics.originalConverged);
    end
end
if ~converged;warning('Cavity:NotConverged','Cavity did not reach its selected stopping criteria.');end

    function [out,distribution,inputWall]=map(state)
        macro=reshape(state(1:16*N),4*N,4);load=zeros(4*N,nv);offset=16*N;
        inputWall.state=cell(1,4);
        for side=1:4
            nodes=op.nodes{side};ids=offset+(1:numel(nodes));offset=offset+numel(nodes);
            incoming=op.cn{side}<0;motion=zeros(1,nv);
            if side==4;motion=2*cfg.lidSpeed*v.vx';end
            incomingState=state(ids)+motion;
            inputWall.state{side}=incomingState;
            load(nodes,incoming)=load(nodes,incoming) ...
                +(op.wallB{side}*incomingState(:,incoming)).*(-op.cn{side}(incoming))';
        end
        rhs=cfg.delta*(op.M*((macro/v.gram)*v.basis'))+load;
        distribution=zeros(4*N,nv);
        for velocity=1:nv;distribution(:,velocity)=op.factor{velocity}\rhs(:,velocity);end
        densityShift=-(op.integral'*(distribution*v.weight))/sum(v.weight);
        distribution=distribution+densityShift;
        for side=1:4;inputWall.state{side}=inputWall.state{side}+densityShift;end
        moments=distribution*v.moment;out=zeros(stateSize,1);out(1:16*N)=moments(:);offset=16*N;
        for side=1:4
            nodes=op.nodes{side};cn=op.cn{side};exiting=cn>0;entering=cn<0;
            density=-distribution(nodes,exiting)*(v.weight(exiting).*cn(exiting))/sum(v.weight(entering).*cn(entering));
            ids=offset+(1:numel(nodes));offset=offset+numel(nodes);out(ids)=density;
        end
        transportCalls=transportCalls+1;
        if mod(transportCalls,cfg.reportEvery)==0
            fprintf('  GMRES transport calls %d, elapsed %.1fs\n',transportCalls,toc(timer));drawnow;
        end
    end

    function [absR,relR,flux,change]=final_diagnostics(distribution,moments,usedWall,gaugeShift)
        current=walls(distribution,mesh,v,op,cfg.lidSpeed);
        gain=(moments/v.gram)*v.basis';residual=-cfg.delta*(op.M*gain)-current.load;
        for velocity=1:nv;residual(:,velocity)=residual(:,velocity)+op.transport{velocity}*distribution(:,velocity);end
        absR=residual_norm(residual,op,v);relR=absR/firstResidual;flux=0;change=0;
        for side=1:4
            cn=op.cn{side};entering=cn<0;exiting=cn>0;nodes=op.nodes{side};u=usedWall.state{side}+gaugeShift;
            mass=distribution(nodes,exiting)*(v.weight(exiting).*cn(exiting))+u(:,entering)*(v.weight(entering).*cn(entering));
            flux=max(flux,max(abs(mass)));change=max(change,max(abs(u(:,entering)-current.state{side}(:,entering)),[],'all'));
        end
    end
end

function n=residual_norm(R,op,v)
n=sqrt(max(0,sum(sum(R.*(op.massFactor\R),1)'.*v.weight)));
end

function density=sweep_wall_density(used,v,op,lidSpeed,gaugeShift)
% Stored incoming amplitudes in the same density gauge as post-sweep H.
density=nan(size(op.M,1),4);
for side=1:4
    rho=used.state{side}(:,1)+gaugeShift;
    if side==4;rho=rho-2*lidSpeed*v.vx(1);end
    density(op.nodes{side},side)=rho;
end
end

function v=velocities(cfg)
for d=1:2
    n=cfg.velocityOrder(d);
    if strcmp(cfg.quadrature,'referenceCubic')
        assert(mod(n,2)==0,'Cubic velocity rule requires even order.');
        z=(2*(1:n)'-1-n)/(n-1);a{d}=cfg.cutoff*z.^3; %#ok<AGROW>
        w{d}=6*cfg.cutoff/(n-1)*z.^2.*exp(-a{d}.^2); %#ok<AGROW>
    else
        b=sqrt((1:n-1)/2);[V,D]=eig(diag(b,1)+diag(b,-1));
        [x,ids]=sort(diag(D));a{d}=(x-flipud(x))/2;
        wt=sqrt(pi)*V(1,ids)'.^2;w{d}=(wt+flipud(wt))/2;
    end
end
[vx,vy]=ndgrid(a{1},a{2});weight=w{1}*w{2}'/pi;
v=struct('vx',vx(:),'vy',vy(:),'weight',weight(:),'normalNodesX',a{1},'normalNodesY',a{2});
v.speed2=v.vx.^2+v.vy.^2;v.basis=[ones(numel(vx),1),v.vx,v.vy,v.speed2-1];
v.moment=v.weight.*v.basis;v.gram=v.basis'*v.moment;
end

function mesh=grid(cells,stretch,mapping)
nx=cells(1);ny=cells(2);
xe=linspace(0,1,nx+1);ye=linspace(0,1,ny+1);
if ismember(mapping,{'quintic','septic','nonic'})
    % Evaluate on the nearer half and reflect, avoiding cancellation at 1.
    sx=min(xe,1-xe);sy=min(ye,1-ye);
    if strcmp(mapping,'quintic')
        mappedX=sx.^3.*(10-15*sx+6*sx.^2);
        mappedY=sy.^3.*(10-15*sy+6*sy.^2);
    elseif strcmp(mapping,'septic')
        mappedX=sx.^4.*(35-84*sx+70*sx.^2-20*sx.^3);
        mappedY=sy.^4.*(35-84*sy+70*sy.^2-20*sy.^3);
    else
        mappedX=sx.^5.*(126-420*sx+540*sx.^2-315*sx.^3+70*sx.^4);
        mappedY=sy.^5.*(126-420*sy+540*sy.^2-315*sy.^3+70*sy.^4);
    end
    mappedX(xe>.5)=1-mappedX(xe>.5);
    mappedY(ye>.5)=1-mappedY(ye>.5);
    xe=mappedX;ye=mappedY;
elseif stretch>0
    xe=.5*(1+tanh(stretch*(2*xe-1))/tanh(stretch));
    ye=.5*(1+tanh(stretch*(2*ye-1))/tanh(stretch));
end
[ix,iy]=ndgrid(1:nx,1:ny);ix=ix(:);iy=iy(:);N=nx*ny;
xl=reshape(xe(ix),[],1);xr=reshape(xe(ix+1),[],1);yb=reshape(ye(iy),[],1);yt=reshape(ye(iy+1),[],1);
dx=xr-xl;dy=yt-yb;
neighbor=[(1:N)'-1,(1:N)'+1,(1:N)'-nx,(1:N)'+nx];
neighbor(ix==1,1)=0;neighbor(ix==nx,2)=0;neighbor(iy==1,3)=0;neighbor(iy==ny,4)=0;
mesh=struct('count',N,'cells',cells,'dx',dx,'dy',dy,'neighbor',neighbor, ...
    'xEdges',xe,'yEdges',ye,'area',dx.*dy, ...
    'center',[(xl+xr)/2,(yb+yt)/2], ...
    'nodes',[reshape([xl,xr,xl,xr]',[],1),reshape([yb,yb,yt,yt]',[],1)]);
end

function op=operators(mesh,v,delta)
N=mesh.count;dx=mesh.dx;dy=mesh.dy;m=[2 1;1 2]/6;g=[-1 -1;1 1]/2;
I=eye(4);E={I(:,[1 3]),I(:,[2 4]),I(:,[1 2]),I(:,[3 4])};op.E=E;
op.M=kron(spdiags(mesh.area,0,N,N),sparse(kron(m,m)));
op.massFactor=decomposition(op.M,'chol');op.integral=repelem(mesh.area/4,4);
Vx=kron(spdiags(dy,0,N,N),sparse(-kron(m,g)));Vy=kron(spdiags(dx,0,N,N),sparse(-kron(g,m)));
op.transport=cell(numel(v.weight),1);op.factor=op.transport;
op.cn={-v.vx,v.vx,-v.vy,v.vy};op.nodes=cell(1,4);op.wallB=cell(1,4);
op.sideLength=[dy dy dx dx];localNodes={[1 3],[2 4],[1 2],[3 4]};op.edgeMass=m;
opposite=[2 1 4 3];C=cell(1,4);B=C;
for s=1:4
    length=op.sideLength(:,s);B{s}=kron(spdiags(length,0,N,N),sparse(E{s}*m*E{s}'));
    ids=find(mesh.neighbor(:,s)>0);neigh=mesh.neighbor(ids,s);
    map=sparse(ids,neigh,length(ids),N,N);
    C{s}=kron(map,sparse(E{s}*m*E{opposite(s)}'));
    boundary=find(mesh.neighbor(:,s)==0);
    op.nodes{s}=reshape((4*(boundary-1)+localNodes{s})',[],1);
    op.wallB{s}=kron(spdiags(length(boundary),0,numel(boundary),numel(boundary)),sparse(m));
end
for j=1:numel(v.weight)
    A=delta*op.M+v.vx(j)*Vx+v.vy(j)*Vy;
    for s=1:4
        A=A+max(op.cn{s}(j),0)*B{s}-max(-op.cn{s}(j),0)*C{s};
    end
    op.transport{j}=A;op.factor{j}=decomposition(A,'lu');
end
end

function wall=walls(H,mesh,v,op,lidSpeed)
wall.load=zeros(size(H));wall.state=cell(1,4);wall.massFlux=0;wall.momentumFlux=zeros(4,2);
for s=1:4
    cn=op.cn{s};incoming=cn<0;outgoing=cn>0;nodes=op.nodes{s};
    denominator=-sum(v.weight(incoming).*cn(incoming));
    motion=zeros(1,numel(v.weight));if s==4;motion=2*lidSpeed*v.vx';end
    rho=-(H(nodes,outgoing)*(v.weight(outgoing).*cn(outgoing)) ...
        +sum(motion(incoming)'.*v.weight(incoming).*cn(incoming)))/sum(v.weight(incoming).*cn(incoming));
    state=rho+motion;wall.state{s}=state;
    wall.load(nodes,incoming)=wall.load(nodes,incoming) ...
        +(op.wallB{s}*state(:,incoming)).*(-cn(incoming))';
    trace=H(nodes,:);trace(:,incoming)=state(:,incoming);
    wall.massFlux=max(wall.massFlux,max(abs(trace*(v.weight.*cn))));
    wall.momentumFlux(s,:)=sum(op.wallB{s},1)*(trace*((v.weight.*cn).*[v.vx,v.vy]));
    assert(denominator>0 && numel(nodes)==2*sum(mesh.neighbor(:,s)==0));
end
end
