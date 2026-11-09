function check_channel_1d2v
%CHECK_CHANNEL_1D2V Independent monolithic DG solve and physical limits.
root = fileparts(fileparts(mfilename('fullpath'))); addpath(root);
diary(fullfile(root,'debug','channel_1d2v_validation.log'));
cleanup = onCleanup(@() diary('off'));
base = struct('plot',false,'save',false,'cells',4,'velocityOrder',[6 4], ...
    'normalQuadrature','hermite','reference',false,'tol',1e-12, ...
    'stopCriterion','equation','reportEvery',1000);
for flow = {'Couette','Poiseuille'}
    opt = base; opt.flow = flow{1};
    r = mono1D2V(opt);
    exact = monolithic_solution(r);
    difference = max(abs(r.H-exact),[],'all');
    assert(r.converged && difference<1e-10);
    assert(isfield(r.mesh,'dx') && ~isfield(r.mesh,'dy'));
    for k = 1:numel(r.velocity.normal)
        assert(all(r.velocity.vx(r.velocity.group{k})==r.velocity.normal(k)));
    end
    assert(max(r.diagnostics.balance)<1e-10);
    assert(r.diagnostics.collisionInvariantDefect<1e-12);
    fprintf('PASS: %s vs independent monolithic DG solve: max H error %.3e.\n',flow{1},difference);
end

% Odd velocity orders include vx=0: these velocities obey collision=source,
% and receive neither an arbitrary transport direction nor a wall inflow.
opt = base; opt.flow = 'Couette'; opt.velocityOrder = [7 5];
opt.wallVelocity = [0.007 0.007]; opt.initial = @(x,vx,vy) 0*x+0*vx+2*0.007*vy;
r = mono1D2V(opt);
assert(r.converged && r.iterations==1);
assert(any(r.velocity.vx==0));
assert(max(abs(r.H-2*0.007*r.velocity.vy'),[],'all')<1e-12);
fprintf('PASS: uniform translating equilibrium and zero normal velocities.\n');

opt = base; opt.cells = 32; opt.velocityOrder = [32 6];
r = mono1D2V(opt);
c = r.couette; p = r.poiseuille;
assert(c.converged && p.converged);
assert(max(abs(c.cellMoments(:,3)+flipud(c.cellMoments(:,3))))<1e-12);
assert(max(abs(p.cellMoments(:,3)-flipud(p.cellMoments(:,3))))<1e-12);
assert(abs(c.flowRate)<1e-12 && p.flowRate>0);
assert(all(diff(c.cellMoments(:,3))>0));
assert(all(p.cellMoments(:,3)>0));
assert(max(abs(c.cellMoments(:,[1 2 4])),[],'all')<1e-12);
assert(max(abs(p.cellMoments(:,[1 2 4])),[],'all')<1e-12);
assert(max(c.diagnostics.balance)<1e-10 && max(p.diagnostics.balance)<1e-10);
assert(abs(c.gasWallVelocity(1))<0.01 && abs(c.gasWallVelocity(2))<0.01);
fprintf('PASS: antisymmetric Couette, symmetric positive Poiseuille, wall slip, conserved face fluxes.\n');

% In the collisionless Couette limit, h=2*vy*U_left for vx>0 and
% h=2*vy*U_right for vx<0. For +/-U walls Pxy=-U/sqrt(pi).
free = base; free.flow = 'Couette'; free.delta = 1e-9;
free.velocityOrder = [32 6];
f = mono1D2V(free);
expected = zeros(size(f.H));
for j = 1:numel(f.velocity.vx)
    side = 1+(f.velocity.vx(j)<0);
    expected(:,j) = 2*f.config.wallVelocity(side)*f.velocity.vy(j);
end
assert(f.converged && max(abs(f.H-expected),[],'all')<1e-9);
discreteStress = expected(1,:)*(f.velocity.weight.*f.velocity.vx.*f.velocity.vy);
assert(max(abs(f.shearStress-discreteStress))<1e-10);
assert(abs(discreteStress+0.01/sqrt(pi))<0.02*0.01/sqrt(pi));
fprintf('PASS: collisionless half-range Couette distribution and shear stress.\n');

% Linear response: reversing the forcing reverses the entire solution.
reverse = opt; reverse.flow = 'Poiseuille'; reverse.acceleration = -0.01;
negative = mono1D2V(reverse);
assert(negative.converged && max(abs(p.H+negative.H),[],'all')<1e-11);
fprintf('PASS: Poiseuille force sign and linear response.\n');

% Ensure tiny moment updates cannot mask a nonzero steady equation residual.
short = base; short.flow = 'Poiseuille'; short.delta = 100; short.maxIter = 1;
short.tol = 1e-4;
old = warning('off','Channel1D2V:NotConverged');
restoreWarning = onCleanup(@() warning(old));
s = mono1D2V(short);
assert(s.history(end,2)<short.tol && ~s.converged && s.diagnostics.kinetic>short.tol);
fprintf('PASS: kinetic residual prevents false convergence from small moment updates.\n');

% Exercise actual plotting and all result files in a separate directory.
restart = p.config; restart.initial = p.H; restart.maxIter = 1;
restart.plot = true; restart.save = true;
restart.outputDir = fullfile(root,'tmp','channel_1d2v_restart');
last = mono1D2V(restart);
folder = fullfile(restart.outputDir,'poiseuille');
assert(last.converged && isfile(fullfile(folder,'summary.png')));
assert(isequal(size(readmatrix(fullfile(folder,'profile.dat'))),[restart.cells 9]));
saved = load(fullfile(folder,'channel.mat'),'r');
assert(saved.r.converged && isequal(last.H,saved.r.H));
close(gcf);
fprintf('PASS: converged restart, plotting, PNG/MAT/DAT/status output.\n');
fprintf('All 1D2V channel checks passed.\n');
end

function exact = monolithic_solution(r)
% Assemble ALL spatial and velocity unknowns together, including current
% collision moments and diffuse-wall density. No lagged production operator.
% Volume matrices are independently integrated by two-point Gauss quadrature.
v = r.velocity; cfg = r.config; N = cfg.cells; n = 2*N; nv = numel(v.vx);
z = [-1 1]/sqrt(3); shape = [(1-z)/2;(1+z)/2];
derivative = [-1;1]/r.mesh.dx;
Mlocal = r.mesh.dx/2*(shape*shape');
Dlocal = r.mesh.dx/2*(derivative*sum(shape,2)');
M = kron(speye(N),sparse(Mlocal));
% Analytic Gaussian invariant norms; independent of production Gram inverse.
P = v.basis*diag([1 2 2 1])*(v.weight.*v.basis)';
A = -cfg.delta*kron(sparse(P),M);
b = zeros(n*nv,1);
for j = 1:nv
    s = v.vx(j); block = (j-1)*n;
    for c = 1:N
        rows = block+2*c-1:block+2*c;
        local = cfg.delta*Mlocal-s*Dlocal;
        if s>0
            local(2,2) = local(2,2)+s;
            if c>1; A(rows(1),rows(1)-1) = A(rows(1),rows(1)-1)-s; end
        elseif s<0
            local(1,1) = local(1,1)-s;
            if c<N; A(rows(2),rows(2)+1) = A(rows(2),rows(2)+1)+s; end
        end
        A(rows,rows) = A(rows,rows)+local;
        b(rows) = Mlocal*ones(2,1)*(2*cfg.acceleration*v.vy(j));
    end
    if s==0; continue; end
    if s>0; node = 1; side = 1; outgoing = find(v.vx<0);
    else; node = n; side = 2; outgoing = find(v.vx>0); end
    row = block+node;
    coeff = v.weight(outgoing).*abs(v.vx(outgoing));
    coeff = coeff/sum(coeff);
    A(row,(outgoing-1)*n+node) = A(row,(outgoing-1)*n+node)-abs(s)*coeff';
    b(row) = b(row)+abs(s)*2*cfg.wallVelocity(side)*v.vy(j);
end
gauge = reshape(ones(n,1)*v.weight'/n,[],1);
full = [A,sparse(gauge);sparse(gauge'),0]\[b;0];
exact = reshape(full(1:end-1),n,nv);
assert(max(abs(A*exact(:)-b))<1e-12);
end
