function check_bgk_dg
%CHECK_BGK_DG Independent DG assembly, wall laws and manufactured solution.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
diary(fullfile(root,'debug','bgk_dg_validation.log'));
cleanup = onCleanup(@() diary('off'));
oldWarning = warning('off','BGKDG:NotConverged');
restoreWarning = onCleanup(@() warning(oldWarning));
base = struct('plot',false,'save',false,'maxIter',1,'velocityOrder',[7 7], ...
    'wallTemperature',[1 1],'initial',0.2,'meanDensity',1.2);
r = monoBGKDG(base);
assert(r.converged && max(abs(r.H-0.2),[],'all')<5e-13);
assert(r.diagnostics.kinetic<1e-12 && r.diagnostics.wall<1e-12);
assert(max(abs(r.mesh.extent-[0 100 0 0.25]))<1e-10);
assert(r.mesh.count==400 && size(r.mesh.faces,1)==1201);
assert(max(abs(r.velocity.gram-diag([1 0.5 0.5 1])),[],'all')<2e-13);
assert(any(r.velocity.vx==0) && any(r.velocity.vy==0));
assert(all(r.velocity.reflect(r.velocity.reflect)==(1:49)'));
fprintf('PASS: original Fluent mesh, odd velocity roots, exact projection, equilibrium and mass constraint.\n');

% Compare one complete production sweep to an independently integrated,
% explicit cell-by-cell implementation of PDF equations (11)-(17).
cfg = struct('plot',false,'save',false,'maxIter',1,'velocityOrder',[7 7], ...
    'wallTemperature',[0.97 1.03],'delta',0.2,'initial',@initial_field);
r = monoBGKDG(cfg);
[reference,stats] = explicit_sweep(r);
difference = max(abs(reference-r.H),[],'all');
assert(difference<5e-12,'Factored sweep differs from direct local DG solves.');
assert(stats.diffuseBalance<1e-13 && stats.crossConstant<1e-13);
assert(r.diagnostics.wallMassFlux>1e-7,'Lagged wall mass-flux monitor must detect intermediate imbalance.');
fprintf('PASS: independent 2x2 Gauss volume/edge assembly and local sweeps, max difference %.3e.\n',difference);
fprintf('PASS: wall quadrature balance %.3e; cross-cell constant trace %.3e; lagged wall flux %.3e.\n', ...
    stats.diffuseBalance,stats.crossConstant,r.diagnostics.wallMassFlux);

% A nontrivial exact piecewise-Q1 solution on the existing 400-cell mesh:
% h*=0.01*tent(x)*(vx^2-1/2)*(vy^2-1/2). All four collision moments vanish,
% it is even in vy, and h*=0 at the two diffuse walls (Tw=1).
manufactured = struct('plot',false,'save',false,'maxIter',1, ...
    'velocityOrder',[8 8],'wallTemperature',[1 1],'delta',0.3, ...
    'initial',@manufactured_h,'source',@manufactured_source);
r = monoBGKDG(manufactured);
xy = r.mesh.points(r.mesh.nodeIds',:);
exact = manufactured_h(xy(:,1),xy(:,2),r.velocity.vx',r.velocity.vy');
difference = max(abs(r.H-exact),[],'all');
assert(r.converged && difference<1e-11 && r.diagnostics.kinetic<1e-11);
assert(r.diagnostics.collisionInvariantDefect<1e-12);
fprintf('PASS: manufactured Q1 steady solution with compatible source, max error %.3e, kinetic %.3e.\n', ...
    difference,r.diagnostics.kinetic);

% Reject a closed-domain source with nonzero total mass production.
bad = manufactured; bad.source = 1;
try
    monoBGKDG(bad);
    error('BGKDG:TestFailure','Incompatible mass source was accepted.');
catch exception
    assert(strcmp(exception.identifier,'BGKDG:SourceMass'));
end
fprintf('PASS: incompatible closed-domain mass source is rejected.\n');

% A stopped iterate must remain visibly unconverged, even if moment changes
% are allowed a loose tolerance. No false success from moments alone.
short = cfg; short.tolMoment = 1;
r = monoBGKDG(short);
assert(~r.converged && r.diagnostics.kinetic>r.config.tolKinetic);
fprintf('PASS: kinetic/wall residuals prevent false convergence from moment changes alone.\n');

folder = fullfile(root,'results','bgk_dg');
if exist(fullfile(folder,'bgk_dg.mat'),'file')
    data = load(fullfile(folder,'bgk_dg.mat'),'result');
    full = data.result;
    assert(full.converged,'Saved full run has not converged.');
    assert(full.diagnostics.kinetic<full.config.tolKinetic);
    assert(full.diagnostics.wall<full.config.tolWall);
    assert(full.diagnostics.wallMassFlux<full.config.tolWall);
    assert(full.diagnostics.massError<full.config.tolWall);
    assert(all(isfinite(full.H),'all'));
    assert(max(abs(full.cellMoments(:,2:3)),[],'all')<1e-7);
    assert(full.cellMoments(1,4)<full.cellMoments(end,4));
    check_saved_fluxes(full);
    % Exercise plotting and MAT/DAT/status output on a converged restart in
    % a separate directory, keeping the original 6989-sweep history intact.
    restart = full.config;
    restart.initial = full.H;
    restart.maxIter = 1;
    restart.plot = true;
    restart.save = true;
    restart.outputDir = fullfile(root,'tmp','bgk_dg_restart');
    last = monoBGKDG(restart);
    assert(last.converged && exist(fullfile(restart.outputDir,'bgk_dg_status.txt'),'file')==2);
    assert(isequal(size(readmatrix(fullfile(restart.outputDir,'bgk_dg_profile.dat'))),[400 8]));
    close(gcf);
    plot_saved_result(full,folder);
    fprintf('PASS: saved two-temperature run, %d sweeps, kinetic %.3e, wall %.3e, mass flux %.3e.\n', ...
        full.iterations,full.diagnostics.kinetic,full.diagnostics.wall,full.diagnostics.wallMassFlux);
else
    fprintf('Full two-temperature convergence run not checked: no saved result yet.\n');
end
fprintf('All BGK DG checks passed.\n');
end

function check_saved_fluxes(r)
% Upwind full-distribution invariant fluxes on INTERNAL vertical faces.
% Their constancy is independent of the solver's residual implementation.
v = r.velocity; H = r.H; N = r.mesh.count;
positive = find(v.vx>0); negative = find(v.vx<0);
right = (H(2:4:end,:)+H(4:4:end,:))/2;
left = (H(1:4:end,:)+H(3:4:end,:))/2;
trace = zeros(N-1,numel(v.vx));
trace(:,positive) = right(1:end-1,positive);
trace(:,negative) = left(2:end,negative);
invariants = [ones(numel(v.vx),1),v.vx,v.vy,v.speed2];
flux = (1+trace)*((v.weight.*v.vx).*invariants);
spread = max(flux,[],1)-min(flux,[],1);
assert(max(spread)<1e-7,'Axial mass/momentum/energy fluxes are not constant.');
fprintf('PASS: axial invariant flux spreads [mass,px,py,energy] = %.3e %.3e %.3e %.3e.\n',spread);
nodal = reshape(H,4,N,[]);
change = nodal(1:2,:,:)-nodal(3:4,:,:);
rawVariation = max(abs(change),[],'all');
% The kinetic stopping norm is M0-weighted. Very small Gaussian tail weights
% do not control the unweighted maximum of h=f/M0-1 at large velocities.
meanSquare = reshape(sum(change.^2,1),N,[]);
dyVariation = sqrt((r.mesh.area'*meanSquare*v.weight)/(2*sum(r.mesh.area)));
assert(dyVariation<1e-8,'A y-independent problem developed a weighted y variation.');
fprintf('PASS: M0-weighted RMS same-x y variation %.3e (unweighted tail maximum %.3e).\n', ...
    dyVariation,rawVariation);
end

function plot_saved_result(r,folder)
f = figure('Visible','off','Color','w','Position',[100 100 1100 650]);
cleanup = onCleanup(@() close(f));
tiledlayout(2,3);
x = r.mesh.center(:,1); q = r.cellMoments;
nexttile; plot(x,q(:,1),'LineWidth',1.3); xlabel('x'); ylabel('density perturbation'); grid on;
nexttile; plot(x,q(:,4),'LineWidth',1.3); xlabel('x'); ylabel('temperature perturbation'); grid on;
nexttile; plot(x,r.heatFlux(:,1),'LineWidth',1.3);
center = mean(r.heatFlux(:,1)); span = max(0.05*abs(center),1e-8);
ylim([center-span,center+span]); xlabel('x'); ylabel('heat flux q_x'); grid on;
nexttile; plot(x,q(:,2:3),'LineWidth',1.2); xlabel('x'); ylabel('velocity moments'); legend('u_x','u_y'); grid on;
nexttile; hold on;
for k = 1:3
    ids = isfinite(r.history(:,k));
    semilogy(r.history(ids,7),max(r.history(ids,k),1e-16),'LineWidth',1.1);
end
set(gca,'YScale','log'); xlabel('source iteration'); ylabel('update / residual'); legend('moments','kinetic','wall'); grid on;
nexttile; hold on;
for k = 4:5
    ids = isfinite(r.history(:,k));
    semilogy(r.history(ids,7),max(r.history(ids,k),eps),'LineWidth',1.1);
end
set(gca,'YScale','log'); xlabel('source iteration'); ylabel('mass flux / mean-density error'); legend('wall flux','mass constraint'); grid on;
exportgraphics(f,fullfile(folder,'bgk_dg_summary.png'),'Resolution',150);
end

function h = initial_field(x,y,vx,vy)
h = 0.02*(x/100+0.2*y/0.25).*(vx.^2-0.5).*(vy.^2-0.5)+0.005*(x/100-0.5);
end

function h = manufactured_h(x,~,vx,vy)
tent = min(x/100,1-x/100);
h = 0.01*tent.*(vx.^2-0.5).*(vy.^2-0.5);
end

function g = manufactured_source(x,y,vx,vy)
center = mean(reshape(x,4,[]),1)';
slope = repelem(sign(50-center)/100,4);
psi = (vx.^2-0.5).*(vy.^2-0.5);
g = 0.01*slope.*vx.*psi+0.3*manufactured_h(x,y,vx,vy);
end

function [next,stats] = explicit_sweep(r)
mesh = r.mesh; vel = r.velocity; cfg = r.config;
N = mesh.count; V = numel(vel.vx);
xy = mesh.points(mesh.nodeIds',:);
H = initial_field(xy(:,1),xy(:,2),vel.vx',vel.vy');
integral = repelem(mesh.area/4,4);
rho = 1+integral'*(H*vel.weight)/sum(mesh.area);
H = (cfg.meanDensity/rho)*(1+H)-1;
moments = H*vel.moment;
gain = (moments/vel.gram)*vel.basis';
[M,Vx,Vy,B,C] = quadrature_matrices(mesh.dx(1),mesh.dy(1));
z = [-1;1]/sqrt(3);
shape = [(1-z)/2,(1+z)/2];
idsL = [1 3]; idsR = 4*(N-1)+[2 4];
left = shape*H(idsL,:); right = shape*H(idsR,:);
rhoL = (1+left(:,vel.outLeft))*(vel.weight(vel.outLeft).*abs(vel.vx(vel.outLeft)))/vel.denominator(1);
rhoR = (1+right(:,vel.outRight))*(vel.weight(vel.outRight).*abs(vel.vx(vel.outRight)))/vel.denominator(2);
gL = rhoL*vel.ratio(:,1)'-1;
gR = rhoR*vel.ratio(:,2)'-1;
left(:,vel.inLeft) = gL(:,vel.inLeft);
right(:,vel.inRight) = gR(:,vel.inRight);
balance = [(1+left)*(vel.weight.*vel.vx);(1+right)*(vel.weight.*vel.vx)];
stats.diffuseBalance = max(abs(balance));
stats.crossConstant = max(abs([C{1}*ones(4,1)-B{1}*ones(4,1);C{2}*ones(4,1)-B{2}*ones(4,1)]));
next = zeros(size(H));
for j = 1:V
    vx = vel.vx(j); vy = vel.vy(j);
    A = cfg.delta*M+vx*Vx+vy*Vy+ ...
        max(vx,0)*B{2}+max(-vx,0)*B{1}+max(vy,0)*B{4}+max(-vy,0)*B{3};
    if vx>=0; sequence = 1:N; else; sequence = N:-1:1; end
    for i = sequence
        ids = 4*(i-1)+(1:4);
        rhs = cfg.delta*M*gain(ids,j);
        if vy>0
            rhs = rhs+vy*B{3}*H(ids,vel.reflect(j));
        elseif vy<0
            rhs = rhs-vy*B{4}*H(ids,vel.reflect(j));
        end
        if vx>0
            if mesh.neighbor(i,1)>0
                up = 4*(mesh.neighbor(i,1)-1)+(1:4);
                rhs = rhs+vx*C{1}*next(up,j);
            else
                rhs([1 3]) = rhs([1 3])+vx*mesh.dy(i)*shape'*gL(:,j)/2;
            end
        elseif vx<0
            if mesh.neighbor(i,2)>0
                up = 4*(mesh.neighbor(i,2)-1)+(1:4);
                rhs = rhs-vx*C{2}*next(up,j);
            else
                rhs([2 4]) = rhs([2 4])-vx*mesh.dy(i)*shape'*gR(:,j)/2;
            end
        end
        next(ids,j) = A\rhs;
    end
end
rho = 1+integral'*(next*vel.weight)/sum(mesh.area);
next = (cfg.meanDensity/rho)*(1+next)-1;
end

function [M,Vx,Vy,B,C] = quadrature_matrices(dx,dy)
% Integrate shape functions directly rather than reuse production matrices.
q = (1+[-1;1]/sqrt(3))/2;
M = zeros(4); Vx = M; Vy = M;
for x = q'
    for y = q'
        phi = [(1-x)*(1-y);x*(1-y);(1-x)*y;x*y];
        gx = [-(1-y);1-y;-y;y]/dx;
        gy = [-(1-x);-x;1-x;x]/dy;
        M = M+(dx*dy/4)*(phi*phi');
        Vx = Vx-(dx*dy/4)*(gx*phi');
        Vy = Vy-(dx*dy/4)*(gy*phi');
    end
end
B = {zeros(4),zeros(4),zeros(4),zeros(4)};
C = {zeros(4),zeros(4)};
for s = q'
    left = [1-s;0;s;0]; right = [0;1-s;0;s];
    bottom = [1-s;s;0;0]; top = [0;0;1-s;s];
    B{1} = B{1}+(dy/2)*(left*left');
    B{2} = B{2}+(dy/2)*(right*right');
    B{3} = B{3}+(dx/2)*(bottom*bottom');
    B{4} = B{4}+(dx/2)*(top*top');
    C{1} = C{1}+(dy/2)*(left*right');
    C{2} = C{2}+(dy/2)*(right*left');
end
end
