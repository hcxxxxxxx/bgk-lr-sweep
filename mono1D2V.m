function result = mono1D2V(options)
%MONO1D2V Standalone linearized BGK / P1-DG channel-flow solver.
%   result = mono1D2V;                % run Couette and Poiseuille separately
%   r = mono1D2V(struct('flow','Couette','delta',10));
%   r = mono1D2V(struct('flow','Poiseuille','acceleration',0.01));
%
%   Physical space: x in [-gap/2,gap/2]; velocities: (vx,vy).
%   x/vx are wall-normal; y/vy are streamwise (no y mesh).
%   M0=exp(-vx^2-vy^2)/pi, f=M0*(1+H), reference rho=T=1.
%   vx*dH/dx = delta*(P(H)-H) + 2*acceleration*vy.
%   P(H)=rho'+2*ux*vx+2*uy*vy+T'*(vx^2+vy^2-1).
%   Both walls are isothermal (T=1), fully diffuse, with tangential speeds
%   wallVelocity=[left,right] along y. Collision gain and wall density lag.
%   A discrete Gram projection conserves all four collision invariants in CIS.
%   iterationMethod='GSIS' uses the reference-style HDG mixed source, without
%   lifting H. monoGSIS1D2V is the separate entry point for that coupling.
%   Default stopping: reference cell-macro relative update < 1e-5, evaluated
%   on mixed source moments in GSIS. Equation residuals remain diagnostics.
%   stopCriterion='equation' restores the earlier residual/update checks.
%   Independent reference: channel1D2VReference.m (disable with reference=false).
%   normalCutoff sets the halfLegendre vx interval [-normalCutoff,normalCutoff].
%   Pure Couette plotting also uses plotCouette1D2V.m for physical checks.
%   No Fluent mesh, monoBGKDG, monoCIS or monoNS dependency.
%
%   Default 'both' returns result.couette and result.poiseuille. A single
%   flow returns its result directly. Run/F5 with no arguments runs both.
%   This is a small-disturbance BGK model, without viscous heating.
if nargin == 0; options = struct; end
cfg = settings(options);
if strcmpi(cfg.flow,'both')
    cfg.flow = 'Couette'; result.couette = solve_channel(case_settings(cfg));
    cfg.flow = 'Poiseuille'; result.poiseuille = solve_channel(case_settings(cfg));
else
    result = solve_channel(case_settings(cfg));
end
end

function cfg = settings(options)
root = fileparts(mfilename('fullpath'));
cfg = struct('flow','both','cells',64,'gap',1,'delta',1, ...
    'iterationMethod','CIS','macroParam',struct, ...
    'velocityOrder',[32 32],'normalQuadrature','halfLegendre','normalCutoff',8, ...
    'wallVelocity',[],'acceleration',[], ...
    'reference',true,'referenceCells',1024,'referenceOrder',160, ...
    'initial',0,'maxIter',30000,'tol',1e-5,'stopCriterion','reference', ...
    'checkEvery',10,'reportEvery',100,'plot',true,'save',true, ...
    'outputDir',fullfile(root,'results','channel_1d2v'));
assert(isstruct(options) && isscalar(options),'options must be a scalar struct.');
names = fieldnames(options);
for k = 1:numel(names)
    assert(isfield(cfg,names{k}),'Unknown option: %s',names{k});
    cfg.(names{k}) = options.(names{k});
end
cfg.flow = char(validatestring(cfg.flow,{'both','Couette','Poiseuille'}));
cfg.iterationMethod = char(validatestring(cfg.iterationMethod,{'CIS','GSIS'}));
cfg.stopCriterion = char(validatestring(cfg.stopCriterion,{'reference','equation'}));
assert(isstruct(cfg.macroParam) && isscalar(cfg.macroParam),'macroParam must be a scalar struct.');
cfg.normalQuadrature = char(validatestring(cfg.normalQuadrature,{'halfLegendre','hermite'}));
for name = {'cells','maxIter','checkEvery','reportEvery','referenceCells','referenceOrder'}
    validateattributes(cfg.(name{1}),{'numeric'},{'scalar','integer','positive','finite'});
end
for name = {'gap','delta','tol','normalCutoff'}
    validateattributes(cfg.(name{1}),{'numeric'},{'scalar','real','positive','finite'});
end
validateattributes(cfg.velocityOrder,{'numeric'},{'vector','numel',2,'integer','>=',3,'finite'});
if strcmp(cfg.normalQuadrature,'halfLegendre')
    assert(mod(cfg.velocityOrder(1),2)==0 && cfg.velocityOrder(1)>=32, ...
        'halfLegendre requires an even Nvx >= 32 (half the points per half-space).');
end
assert(cfg.referenceOrder>=4,'referenceOrder must be >= 4.');
if ~isempty(cfg.wallVelocity)
    validateattributes(cfg.wallVelocity,{'numeric'},{'vector','numel',2,'real','finite'});
    cfg.wallVelocity = cfg.wallVelocity(:)';
end
if ~isempty(cfg.acceleration)
    validateattributes(cfg.acceleration,{'numeric'},{'scalar','real','finite'});
end
for name = {'plot','save','reference'}
    validateattributes(cfg.(name{1}),{'logical','numeric'},{'scalar','binary'});
end
end

function cfg = case_settings(cfg)
if isempty(cfg.wallVelocity)
    if strcmpi(cfg.flow,'Couette'); cfg.wallVelocity = [-0.01 0.01];
    else; cfg.wallVelocity = [0 0]; end
end
if isempty(cfg.acceleration)
    if strcmpi(cfg.flow,'Poiseuille'); cfg.acceleration = 0.01;
    else; cfg.acceleration = 0; end
end
end

function r = solve_channel(cfg)
N = cfg.cells; dx = cfg.gap/N;
edges = linspace(-cfg.gap/2,cfg.gap/2,N+1)';
x = reshape([edges(1:end-1),edges(2:end)]',[],1);
mesh = struct('count',N,'edges',edges,'nodes',x, ...
    'center',(edges(1:end-1)+edges(2:end))/2,'dx',dx);
v = velocity_grid(cfg.velocityOrder,cfg.normalQuadrature,cfg.normalCutoff);
op = operators(N,dx,v,cfg.delta);
ref = [];
if cfg.reference
    ref = channel1D2VReference(cfg,ceil(cfg.referenceCells/N)*N,cfg.referenceOrder);
    ref.cellVelocity = mean(reshape(ref.velocity,ref.cells/N,N),1)';
end
if isa(cfg.initial,'function_handle')
    H = cfg.initial(x,v.vx',v.vy');
else
    H = cfg.initial;
end
if isscalar(H); H = repmat(H,2*N,numel(v.vx)); end
assert(isequal(size(H),[2*N,numel(v.vx)]) && isreal(H) && all(isfinite(H),'all'), ...
    'initial must be a scalar, 2*Ncell by Nvelocity array, or @(x,vx,vy).');
% Additive density gauge: the domain-average density perturbation is zero.
H = H-mean(H*v.weight);
G = repmat(2*cfg.acceleration*v.vy',2*N,1);
Q = H*v.moment;
sourceQ=Q; useGSIS=strcmpi(cfg.iterationMethod,'GSIS');
referenceColumns=7+2*useGSIS+(1:4);
history = nan(cfg.maxIter,referenceColumns(end)); converged = false; monitor = [];
previousReference=[diff(edges),diff(edges),zeros(N,1)];
macroParam=cfg.macroParam; synthetic=struct;
if useGSIS
    macroParam.cache=hdg_precompute_p1(mesh,cfg.delta,macroParam);
    omega=repelem(min(1,2*cfg.delta*diff(edges)/sqrt(pi)),2);
end
if cfg.plot
    monitor = figure('Name',['1D2V ',cfg.flow],'Color','w','Position',[100 100 1100 720]);
    if isprop(monitor,'Theme'); monitor.Theme = 'light'; end
end
fprintf('\n1D2V %s %s: P1-DG, %d cells, %dx%d velocities, gap=%g, delta=%g\n', ...
    cfg.iterationMethod,cfg.flow,N,cfg.velocityOrder,cfg.gap,cfg.delta);
fprintf('x-normal mesh; wall uy=[%g %g], y-acceleration=%g; fully diffuse walls.\n', ...
    cfg.wallVelocity,cfg.acceleration);
fprintf('Normal quadrature: %s.\n',cfg.normalQuadrature);
fprintf('Stopping: %s, tolerance %.3e.\n',cfg.stopCriterion,cfg.tol);
fprintf(' step       moment update     kinetic residual    wall residual      wall mass flux\n');
if useGSIS
    fprintf('GSIS S=[%s], omega=[%.4g,%.4g]; additional column: coupled residual\n', ...
        num2str(macroParam.cache.stabilization),min(omega),max(omega));
end
timer = tic;
for iter = 1:cfg.maxIter
    oldQ = Q;
    oldSourceQ=sourceQ;
    used = wall_data(H,v,cfg.wallVelocity);
    gain = (sourceQ/v.gram)*v.basis';
    rhs = op.M*(cfg.delta*gain+G)+wall_load(used,v,2*N);
    % Each vx shares one sparse, block-triangular transport system for all vy.
    for k = 1:numel(v.normal)
        ids = v.group{k};
        H(:,ids) = op.factor{k}\rhs(:,ids);
    end
    shift = -mean(H*v.weight);
    H = H+shift; used = used+shift;
    Q = H*v.moment;
    if useGSIS
        state=struct('mesh',mesh,'H',H,'velocity',v,'config',cfg);
        [macro,aux,hdg,hot]=gsis_macro_from_kinetic_1d2v(state,macroParam);
        sourceQ=Q+omega.*(macro.nodal-Q);
        sourceQ(:,1)=sourceQ(:,1)-mean(sourceQ(:,1));
        synthetic=struct('macro',macro,'aux',aux,'hdg',hdg,'hot',hot,'omega',omega);
    else
        sourceQ=Q;
    end
    sourceUpdate=max(sqrt(mean((sourceQ-oldSourceQ).^2,1)));
    update = max(sqrt(mean((Q-oldQ).^2,1)));
    [referenceUpdate,referenceComponents,previousReference]= ...
        reference_macro_update_1d2v(sourceQ,edges,previousReference);
    history(iter,referenceColumns)=[referenceUpdate,referenceComponents];
    check = iter==1 || mod(iter,cfg.checkEvery)==0 || ...
        mod(iter,cfg.reportEvery)==0 || iter==cfg.maxIter || update<cfg.tol ...
        || referenceUpdate<cfg.tol;
    if check
        d = diagnostics(H,Q,G,used,v,op,cfg,mesh,sourceQ);
        if strcmpi(cfg.stopCriterion,'reference')
            converged=referenceUpdate<cfg.tol;
        else
            converged=max([update,sourceUpdate,d.coupled,d.wall,d.wallMassFlux,d.massError])<cfg.tol;
        end
        history(iter,1:7) = [iter,update,d.kinetic,d.wall,d.wallMassFlux,d.massError,toc(timer)];
        if useGSIS; history(iter,8:9)=[sourceUpdate,d.coupled]; end
    else
        history(iter,[1 2 7]) = [iter,update,toc(timer)];
    end
    if iter==1 || mod(iter,cfg.reportEvery)==0 || converged || iter==cfg.maxIter
        fprintf('%6d      %.5e       %.5e       %.5e       %.5e\n',iter,history(iter,2:5));
        if useGSIS; fprintf('            source update %.5e; coupled residual %.5e\n',sourceUpdate,d.coupled); end
        fprintf('            reference macro update %.5e [rho %.3e, T %.3e, speed %.3e]\n', ...
            referenceUpdate,referenceComponents);
        if isgraphics(monitor); plot_progress(monitor,mesh,Q,history(1:iter,:),cfg,H,v,ref); end
        drawnow;
    end
    if converged; break; end
end
r = struct('config',cfg,'mesh',mesh,'velocity',v,'H',H,'nodalMoments',Q, ...
    'cellMoments',cell_average(Q),'history',history(1:iter,:), ...
    'converged',converged,'iterations',iter,'wallTime',toc(timer));
r.historyColumns = {'iteration','momentUpdate','kineticResidual', ...
    'wallResidual','wallMassFlux','massError','elapsedSeconds'};
if useGSIS
    r.historyColumns=[r.historyColumns,{'sourceUpdate','coupledResidual'}];
    synthetic.sourceMoments=sourceQ;
    r.synthetic=synthetic;
end
r.historyColumns=[r.historyColumns,{'referenceMacroUpdate','relativeDensityUpdate', ...
    'relativeTemperatureUpdate','relativeSpeedUpdate'}];
r.convergence=struct('criterion',cfg.stopCriterion,'tolerance',cfg.tol, ...
    'referenceMacroUpdate',referenceUpdate,'relativeMacroComponents',referenceComponents, ...
    'cellMonitor',previousReference);
r.momentColumns = {'rhoPerturbation','ux','uy','temperaturePerturbation'};
r.coordinates = struct('space','x','normalVelocity','vx','streamwiseVelocity','vy');
r.shearStress = cell_average(H*(v.weight.*v.vx.*v.vy));
r.heatFlux = cell_average(H*(v.weight.*[v.vx.*(v.speed2-2),v.vy.*(v.speed2-2)]));
r.kineticTrace = face_traces(H,used,v,mesh.count);
r.flowRate = dx*sum(r.cellMoments(:,3));
r.noSlipVelocity = no_slip(mesh.center,cfg);
r.diagnostics = diagnostics(H,Q,G,used,v,op,cfg,mesh,sourceQ);
% Gas-side endpoint moments (DG traces); not prescribed wall speeds.
r.gasWallVelocity = Q([1 end],3)';
r.slipVelocity = r.gasWallVelocity-cfg.wallVelocity;
r.reference = ref;
if cfg.reference
    scale = norm(ref.cellVelocity);
    if scale>0; r.reference.relativeL2 = norm(r.cellMoments(:,3)-ref.cellVelocity)/scale;
    else; r.reference.relativeL2 = nan; end
    r.reference.maximumAbsoluteError = max(abs(r.cellMoments(:,3)-ref.cellVelocity));
end
if cfg.plot && isgraphics(monitor) && strcmpi(cfg.flow,'Couette') && cfg.acceleration==0 && N>=4
    r.couetteChecks = plotCouette1D2V(r,cfg.save);
end
if cfg.save; write_results(r,monitor); end
fprintf('%s finished: converged=%d, %d sweeps, %.2fs, flow rate %.9g\n', ...
    cfg.flow,converged,iter,r.wallTime,r.flowRate);
if cfg.reference
    fprintf('Same-model integral reference: velocity relative L2 %.5e, reference flow rate %.10g\n', ...
        r.reference.relativeL2,ref.flowRate);
end
if ~converged
    warning('Channel1D2V:NotConverged', ...
        '%s %s reached maxIter; reference update %.3e, coupled residual %.3e.', ...
        cfg.iterationMethod,cfg.flow,referenceUpdate,r.diagnostics.coupled);
end
end

function v = velocity_grid(order,normalQuadrature,normalCutoff)
if strcmp(normalQuadrature,'hermite')
    [x,wx] = hermite_rule(order(1));
else
    % Separate half spaces resolve the steep behavior near vx=0 at walls.
    % Include M0's normal Gaussian in the weights; default cutoff is 8.
    n = order(1)/2; k = (1:n-1)'; b = k./sqrt(4*k.^2-1);
    [V,D] = eig(diag(b,1)+diag(b,-1));
    [z,idx] = sort(diag(D)); s = (normalCutoff/2)*(z+1);
    w = normalCutoff*V(1,idx)'.^2.*exp(-s.^2);
    x = [-flipud(s);s]; wx = [flipud(w);w];
end
[z,wz] = hermite_rule(order(2));
[xx,zz] = meshgrid(x,z); w = wz*wx'/pi;
v.vx = xx(:); v.vy = zz(:); v.weight = w(:);
v.speed2 = v.vx.^2+v.vy.^2;
v.basis = [ones(numel(xx),1),v.vx,v.vy,v.speed2-1];
v.moment = v.weight.*v.basis; v.gram = v.basis'*v.moment;
v.normal = x; v.group = cell(numel(x),1);
for k = 1:numel(x); v.group{k} = (k-1)*numel(z)+(1:numel(z)); end
v.inLeft = v.vx>0; v.inRight = v.vx<0;
v.halfFlux = [sum(v.weight.*max(v.vx,0)),sum(v.weight.*max(-v.vx,0))];
v.normalQuadrature = normalQuadrature;
v.normalCutoff = normalCutoff;
end

function [x,w] = hermite_rule(n)
a = sqrt((1:n-1)/2);
[U,D] = eig(diag(a,1)+diag(a,-1));
[x,idx] = sort(diag(D)); w = sqrt(pi)*U(1,idx)'.^2;
x = (x-flipud(x))/2; w = (w+flipud(w))/2;
x(abs(x)<1e-14) = 0;
end

function op = operators(N,dx,v,delta)
m = dx*[2 1;1 2]/6;
d = [-1 -1;1 1]/2; % integral (d test/dx)*trial
op.M = kron(speye(N),sparse(m));
op.massFactor = decomposition(op.M,'chol');
op.transport = cell(numel(v.normal),1); op.factor = op.transport;
for k = 1:numel(v.normal)
    s = v.normal(k);
    local = delta*m-s*d+diag([max(-s,0),max(s,0)]);
    A = kron(speye(N),sparse(local));
    if s>0
        A = A+sparse(3:2:2*N,2:2:2*N-2,-s*ones(1,N-1),2*N,2*N);
    elseif s<0
        A = A+sparse(2:2:2*N-2,3:2:2*N,s*ones(1,N-1),2*N,2*N);
    end
    op.transport{k} = A; op.factor{k} = decomposition(A,'lu');
end
end

function wall = wall_data(H,v,U)
% h_in = rho_wall' + 2*vy*Uwall; rho_wall' enforces zero x mass flux.
leftOut = v.vx<0; rightOut = v.vx>0;
rhoL = H(1,leftOut)*(v.weight(leftOut).*abs(v.vx(leftOut)))/v.halfFlux(1);
rhoR = H(end,rightOut)*(v.weight(rightOut).*abs(v.vx(rightOut)))/v.halfFlux(2);
wall = [rhoL+2*U(1)*v.vy';rhoR+2*U(2)*v.vy'];
end

function load = wall_load(wall,v,n)
load = zeros(n,numel(v.vx));
load(1,v.inLeft) = wall(1,v.inLeft).*v.vx(v.inLeft)';
load(end,v.inRight) = -wall(2,v.inRight).*v.vx(v.inRight)';
end

function d = diagnostics(H,Q,G,used,v,op,cfg,mesh,sourceQ)
wall = wall_data(H,v,cfg.wallVelocity);
gain = (Q/v.gram)*v.basis';
R = -op.M*(cfg.delta*gain+G)-wall_load(wall,v,size(H,1));
for k = 1:numel(v.normal)
    ids = v.group{k}; R(:,ids) = R(:,ids)+op.transport{k}*H(:,ids);
end
d.kinetic = sqrt(max(0,sum(sum(R.*(op.massFactor\R),1)'.*v.weight)/cfg.gap));
Rsource=R-op.M*(cfg.delta*((sourceQ-Q)/v.gram)*v.basis');
d.coupled=sqrt(max(0,sum(sum(Rsource.*(op.massFactor\Rsource),1)'.*v.weight)/cfg.gap));
d.sourceMomentGap=max(sqrt(mean((sourceQ-Q).^2,1)));
change = (wall-used).^2;
d.wall = sqrt(max([change(1,v.inLeft)*v.weight(v.inLeft), ...
    change(2,v.inRight)*v.weight(v.inRight)]));
trace = face_traces(H,used,v,mesh.count);
% Physical numerical face fluxes, with the ACTUAL last-sweep wall inflows.
d.faceFlux = (1+trace)*((v.weight.*v.vx).* ...
    [ones(numel(v.vx),1),v.vy,v.vx,v.speed2]);
d.faceFluxColumns = {'mass','streamwiseMomentum','normalMomentum','energy'};
d.wallMassFlux = max(abs(d.faceFlux([1 end],1)));
d.massError = abs(mean(Q(:,1)));
d.collisionInvariantDefect = max(abs((H-gain)*v.moment),[],'all');
% d(Pxy)/dx=a. All other invariant fluxes are constant in these tests.
balance = diff(d.faceFlux,1,1)/mesh.dx-repmat([0 cfg.acceleration 0 0],mesh.count,1);
d.balance = max(abs(balance),[],1);
end

function trace = face_traces(H,wall,v,N)
trace = zeros(N+1,numel(v.vx));
trace(1,:) = H(1,:); trace(end,:) = H(end,:);
trace(1,v.inLeft) = wall(1,v.inLeft);
trace(end,v.inRight) = wall(2,v.inRight);
trace(2:N,v.inLeft) = H(2:2:end-2,v.inLeft);
trace(2:N,v.inRight) = H(3:2:end,v.inRight);
end

function q = cell_average(Q)
q = (Q(1:2:end,:)+Q(2:2:end,:))/2;
end

function u = no_slip(x,cfg)
z = (x+cfg.gap/2)/cfg.gap;
% M0 has p=1/2 and mu=p/delta=1/(2*delta).
u = cfg.wallVelocity(1)*(1-z)+cfg.wallVelocity(2)*z + ...
    cfg.acceleration*cfg.delta*(x+cfg.gap/2).*(cfg.gap/2-x);
end

function plot_progress(window,mesh,Q,history,cfg,H,v,ref)
figure(window); clf(window); tiledlayout(2,2);
q = cell_average(Q);
nexttile; plot(mesh.center,q(:,3),'LineWidth',1.3); hold on;
labels = {'P1-DG'};
if ~isempty(ref)
    plot(mesh.center,ref.cellVelocity,'--','LineWidth',1.3);
    labels{end+1} = 'BGK integral reference';
end
plot(mesh.center,no_slip(mesh.center,cfg),':','Color',[0.5 0.5 0.5],'LineWidth',1.1);
plot([-cfg.gap/2 cfg.gap/2],cfg.wallVelocity,'ks');
labels = [labels,{'NS no slip (continuum limit)','walls'}];
xlabel('x'); ylabel('u_y'); legend(labels,'Location','best'); grid on;
nexttile;
if ~isempty(ref)
    plot(mesh.center,q(:,3)-ref.cellVelocity,'LineWidth',1.2);
    xlabel('x'); ylabel('u_y(DG) - u_y(BGK reference)'); grid on;
else
    plot(mesh.center,q(:,[1 2 4]),'LineWidth',1.1);
    xlabel('x'); legend('rho perturbation','u_x','T perturbation','Location','best'); grid on;
    span = max(1e-12,1.1*max(abs(q(:,[1 2 4])),[],'all')); ylim([-span span]);
end
nexttile; hold on;
columns=2:4; labelsResidual={'moments','kinetic CIS','wall'};
if strcmpi(cfg.iterationMethod,'GSIS')
    columns=[columns 8 9]; labelsResidual=[labelsResidual,{'source update','coupled'}];
end
columns=[columns,size(history,2)-3]; labelsResidual=[labelsResidual,{'reference macro update'}];
for k = columns
    ids = isfinite(history(:,k));
    semilogy(history(ids,1),max(history(ids,k),1e-16),'.-','LineWidth',1.1);
end
set(gca,'YScale','log'); xlabel('source iteration'); ylabel('update / residual');
legend(labelsResidual); grid on;
nexttile;
stress = cell_average(H*(v.weight.*v.vx.*v.vy));
plot(mesh.center,stress,'LineWidth',1.2); hold on;
plot(mesh.center,mean(stress)+cfg.acceleration*mesh.center,'--','LineWidth',1.1);
xlabel('x'); ylabel('P_{xy}'); legend('BGK','mean(P_{xy}) + a x','Location','best'); grid on;
span = max([0.05*abs(mean(stress)),0.6*(max(stress)-min(stress)),1e-8]);
ylim(mean(stress)+[-span span]);
if ~isempty(ref) && norm(ref.cellVelocity)>0
    err = norm(q(:,3)-ref.cellVelocity)/norm(ref.cellVelocity);
    sgtitle(sprintf('1D2V %s, delta=%g; velocity relative L2 vs BGK reference = %.2e',cfg.flow,cfg.delta,err));
else
    sgtitle(sprintf('1D2V %s, delta=%g',cfg.flow,cfg.delta));
end
drawnow;
end

function write_results(r,monitor)
folder = fullfile(r.config.outputDir,lower(r.config.flow));
if ~exist(folder,'dir'); mkdir(folder); end
save(fullfile(folder,'channel.mat'),'r');
writematrix([r.mesh.center,r.cellMoments,r.shearStress,r.heatFlux,r.noSlipVelocity], ...
    fullfile(folder,'profile.dat'),'Delimiter','tab');
writematrix(r.history,fullfile(folder,'history.dat'),'Delimiter','tab');
writematrix([r.mesh.edges,r.diagnostics.faceFlux],fullfile(folder,'face_flux.dat'),'Delimiter','tab');
fid = fopen(fullfile(folder,'status.txt'),'w'); assert(fid>=0,'Cannot write status file.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'converged=%d\niterations=%d\nwall_time=%.9g\nflow_rate=%.12e\n', ...
    r.converged,r.iterations,r.wallTime,r.flowRate);
fprintf(fid,'coordinates: mesh=x normal=vx streamwise=vy; wallVelocity=[left right] along y\n');
fprintf(fid,'kinetic_residual=%.12e\nwall_residual=%.12e\nwall_mass_flux=%.12e\n', ...
    r.diagnostics.kinetic,r.diagnostics.wall,r.diagnostics.wallMassFlux);
fprintf(fid,'iteration_method=%s\ncoupled_residual=%.12e\nsource_moment_gap=%.12e\n', ...
    r.config.iterationMethod,r.diagnostics.coupled,r.diagnostics.sourceMomentGap);
fprintf(fid,'stop_criterion=%s\nstop_tolerance=%.12e\nreference_macro_update=%.12e\n', ...
    r.convergence.criterion,r.convergence.tolerance,r.convergence.referenceMacroUpdate);
fprintf(fid,'relative_density_update=%.12e\nrelative_temperature_update=%.12e\nrelative_speed_update=%.12e\n', ...
    r.convergence.relativeMacroComponents);
fprintf(fid,'gas_wall_uy=%.12e %.12e\nslip_uy=%.12e %.12e\n',r.gasWallVelocity,r.slipVelocity);
if ~isempty(r.reference)
    writematrix([r.mesh.center,r.cellMoments(:,3),r.reference.cellVelocity, ...
        r.cellMoments(:,3)-r.reference.cellVelocity],fullfile(folder,'reference_comparison.dat'),'Delimiter','tab');
    fprintf(fid,'BGK_reference_relative_L2=%.12e\nBGK_reference_max_absolute_error=%.12e\n', ...
        r.reference.relativeL2,r.reference.maximumAbsoluteError);
    fprintf(fid,'reference_comparison: x uy_DG uy_BGK_reference difference\n');
end
fprintf(fid,'profile: x rho_pert ux uy T_pert Pxy heat_x heat_y uy_NS_no_slip\n');
fprintf(fid,'history: %s\n',strjoin(r.historyColumns,' '));
fprintf(fid,'face_flux: x mass streamwise_momentum normal_momentum energy\n');
if isgraphics(monitor); exportgraphics(monitor,fullfile(folder,'summary.png'),'Resolution',150); end
end
