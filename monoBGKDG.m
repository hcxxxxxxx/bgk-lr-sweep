function result = monoBGKDG(options)
%MONOBGKDG 2D2V linearized BGK, Q1 upwind DG and lagged source iteration.
%   result = monoBGKDG;                       % existing mesh/nsw400.cas
%   result = monoBGKDG(struct('delta',0.1,'wallTemperature',[0.98 1.02]));
%   Node order: bottom-left, bottom-right, top-left, top-right.
%   H is (4*Ncell)-by-Nvelocity; moments are perturbations, NOT NS states.
%   See README.md, section "PDF BGK / Q1-DG solver" for options and output.
if nargin == 0; options = struct; end
cfg = settings(options);
mesh = read_strip(cfg);
vel = velocity_grid(cfg.velocityOrder,cfg.wallTemperature,cfg.wallModel);
op = operators(mesh,vel,cfg.delta);
[x,y] = nodal_coordinates(mesh);
H = field_values(cfg.initial,x,y,vel);
G = field_values(cfg.source,x,y,vel);
sourceLoad = op.M*G;
volume = sum(mesh.area);
sourceMass = op.integral'*(G*vel.weight);
if abs(sourceMass) > 1e-11*max(1,volume*max(abs(G),[],'all'))
    error('BGKDG:SourceMass','A closed steady domain requires integral(G*M0)=0.');
end
sourceFree = ~any(G,'all');
H = constrain_mass(H,op,vel,cfg.meanDensity,sourceFree);
Q = H*vel.moment;
history = nan(cfg.maxIter,7);
converged = false;
timer = tic;
fprintf('BGK Q1-DG: %d cells, %d velocities, domain %.6g x %.6g, delta %.6g\n', ...
    mesh.count,numel(vel.vx),mesh.extent(2)-mesh.extent(1), ...
    mesh.extent(4)-mesh.extent(3),cfg.delta);
fprintf('Tw=[%.6g %.6g] (%s wall); incoming physical-wall data and collision gain are lagged.\n', ...
    cfg.wallTemperature,cfg.wallModel);
fprintf(' step       moment update      kinetic residual      wall residual       wall mass flux\n');
monitor = [];
if cfg.plot
    monitor = figure('Name','Linearized BGK Q1-DG','Color','w');
end
for iter = 1:cfg.maxIter
    previousQ = Q;
    % All velocity moments and physical-wall traces belong to iteration n.
    gain = (Q/vel.gram)*vel.basis';
    wall = wall_data(H,mesh,vel,op);
    rhs = cfg.delta*(op.M*gain) + sourceLoad + wall.load;
    for j = 1:numel(vel.vx)
        % This is a factored block-triangular upwind sweep, with current
        % upstream cell traces. Each velocity system has four unknowns/cell.
        H(:,j) = op.factor{j}\rhs(:,j);
    end
    [H,affine] = constrain_mass(H,op,vel,cfg.meanDensity,sourceFree);
    usedWall = struct('gLeft',affine(1)*wall.gLeft+affine(2), ...
        'gRight',affine(1)*wall.gRight+affine(2), ...
        'reflectedH',affine(1)*wall.reflectedH+affine(2));
    Q = H*vel.moment;
    update = max(sqrt(sum((Q-previousQ).^2,1)./max(sum(previousQ.^2,1),size(Q,1))));
    check = iter == 1 || mod(iter,cfg.checkEvery) == 0 || ...
        mod(iter,cfg.reportEvery) == 0 || iter == cfg.maxIter || update < cfg.tolMoment;
    kinetic = nan; reflection = nan; massFlux = nan; massError = nan;
    if check
        diagnostics = steady_diagnostics(H,Q,G,mesh,vel,op,cfg,usedWall);
        kinetic = diagnostics.kinetic;
        reflection = diagnostics.wall;
        massFlux = diagnostics.wallMassFlux;
        massError = diagnostics.massError;
        converged = update < cfg.tolMoment && kinetic < cfg.tolKinetic && ...
            reflection < cfg.tolWall && massFlux < cfg.tolWall && massError < cfg.tolWall;
    end
    history(iter,:) = [update,kinetic,reflection,massFlux,massError,toc(timer),iter];
    if iter == 1 || mod(iter,cfg.reportEvery) == 0 || converged || iter == cfg.maxIter
        fprintf('%6d      %.6e       %.6e        %.6e       %.6e\n', ...
            iter,update,kinetic,reflection,massFlux);
        if cfg.plot && isgraphics(monitor)
            plot_progress(monitor,mesh,Q,history(1:iter,:));
        end
        drawnow;
    end
    if converged; break; end
end
diagnostics = steady_diagnostics(H,Q,G,mesh,vel,op,cfg,usedWall);
result.config = cfg;
result.mesh = mesh;
result.velocity = vel;
result.H = H;
result.moments = reshape(Q,4,mesh.count,4);
result.cellMoments = reshape(mean(result.moments,1),mesh.count,4);
result.history = history(1:iter,:);
result.historyColumns = {'momentUpdate','kineticResidual','wallResidual', ...
    'wallMassFlux','massError','elapsedSeconds','iteration'};
result.diagnostics = diagnostics;
result.converged = converged;
result.iterations = iter;
result.wallTime = toc(timer);
result.heatFlux = cell_average(H*(vel.weight.* ...
    [vel.vx.*(vel.speed2-2),vel.vy.*(vel.speed2-2)]),mesh.count);
result.minimumFullDistributionRatio = min(1+H,[],'all');
if cfg.save
    write_results(result);
end
if converged
    fprintf('BGK DG converged in %d sweeps (%.2fs), kinetic %.3e, wall %.3e.\n', ...
        iter,result.wallTime,diagnostics.kinetic,diagnostics.wall);
else
    warning('BGKDG:NotConverged', ...
        'Stopped at %d sweeps; kinetic %.3e, wall %.3e. result.converged=false.', ...
        iter,diagnostics.kinetic,diagnostics.wall);
end
end

function cfg = settings(options)
root = fileparts(mfilename('fullpath'));
cfg = struct('meshFile',fullfile(root,'mesh','nsw400.cas'), ...
    'domain',[],'delta',0.01,'wallTemperature',[0.95 1.05], ...
    'wallModel','nonlinear','velocityOrder',[16 16],'meanDensity',1, ...
    'source',0,'initial',0, ...
    'maxIter',30000,'tolMoment',1e-9,'tolKinetic',1e-9,'tolWall',1e-9, ...
    'checkEvery',20,'reportEvery',200,'plot',true,'save',true, ...
    'outputDir',fullfile(root,'results','bgk_dg'));
assert(isstruct(options) && isscalar(options),'options must be a scalar struct.');
names = fieldnames(options);
for k = 1:numel(names)
    assert(isfield(cfg,names{k}),'Unknown option: %s',names{k});
    cfg.(names{k}) = options.(names{k});
end
validateattributes(cfg.delta,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(cfg.wallTemperature,{'numeric'},{'vector','numel',2,'real','finite','positive'});
cfg.wallModel = char(validatestring(cfg.wallModel,{'nonlinear','linear'}));
validateattributes(cfg.velocityOrder,{'numeric'},{'vector','numel',2,'integer','>=',3});
validateattributes(cfg.meanDensity,{'numeric'},{'scalar','real','finite','positive'});
for name = {'maxIter','checkEvery','reportEvery'}
    validateattributes(cfg.(name{1}),{'numeric'},{'scalar','integer','positive','finite'});
end
for name = {'tolMoment','tolKinetic','tolWall'}
    validateattributes(cfg.(name{1}),{'numeric'},{'scalar','positive','finite'});
end
if ~isempty(cfg.domain)
    validateattributes(cfg.domain,{'numeric'},{'vector','numel',4,'real','finite'});
    assert(cfg.domain(2)>cfg.domain(1) && cfg.domain(4)>cfg.domain(3));
end
end

function mesh = read_strip(cfg)
% Recover every rectangle and neighbor from the supplied Fluent topology.
b = importFluentMesh2(cfg.meshFile,1);
p = b.points;
if ~isempty(cfg.domain)
    for d = 1:2
        old = [min(p(:,d)),max(p(:,d))];
        target = cfg.domain(2*d-1:2*d);
        p(:,d) = target(1)+(p(:,d)-old(1))*diff(target)/diff(old);
    end
end
N = b.meshInfo(3);
boxes = zeros(N,4);
nodes = cell(N,1);
cellFaces = cell(N,1);
for e = 1:size(b.faces,1)
    for c = b.faces(e,3:4)
        if c > 0; cellFaces{c}(end+1) = e; end
    end
end
tol = 1e-8*max(1,max(p,[],'all')-min(p,[],'all'));
for c = 1:N
    ids = unique(reshape(b.faces(cellFaces{c},1:2),[],1));
    xy = p(ids,:);
    assert(numel(ids)==4 && numel(cellFaces{c})==4, ...
        'BGKDG requires one row of axis-aligned rectangular cells.');
    boxes(c,:) = [min(xy(:,1)),max(xy(:,1)),min(xy(:,2)),max(xy(:,2))];
    box = boxes(c,:);
    corners = [box(1) box(3);box(2) box(3);box(1) box(4);box(2) box(4)];
    node = zeros(1,4);
    for a = 1:4
        hit = find(max(abs(xy-corners(a,:)),[],2)<tol);
        assert(isscalar(hit),'Cell %d is not an axis-aligned rectangle.',c);
        node(a) = ids(hit);
    end
    nodes{c} = node;
end
[~,order] = sort(boxes(:,1));
inverse(order) = 1:N;
boxes = boxes(order,:);
extent = [min(p(:,1)),max(p(:,1)),min(p(:,2)),max(p(:,2))];
assert(all(abs(boxes(:,3)-extent(3))<tol) && all(abs(boxes(:,4)-extent(4))<tol), ...
    'The PDF solver requires Ny=1.');
assert(all(abs(boxes(1:end-1,2)-boxes(2:end,1))<tol),'The strip has gaps or overlaps.');
neighbor = zeros(N,4); faceIds = zeros(N,4);
for i = 1:N
    original = order(i);
    box = boxes(i,:);
    center = [(box(1)+box(2))/2,(box(3)+box(4))/2];
    targets = [box(1) center(2);box(2) center(2);center(1) box(3);center(1) box(4)];
    for e = cellFaces{original}
        midpoint = mean(p(b.faces(e,1:2),:),1);
        side = find(max(abs(targets-midpoint),[],2)<tol);
        assert(isscalar(side),'Invalid face in rectangle %d.',original);
        faceIds(i,side) = e;
        adjacent = b.faces(e,3:4);
        other = adjacent(adjacent~=original & adjacent>0);
        if ~isempty(other); neighbor(i,side) = inverse(other); end
    end
end
assert(all(neighbor(:,3:4)==0,'all') && neighbor(1,1)==0 && neighbor(end,2)==0);
assert(isequal(neighbor(2:end,1),(1:N-1)') && isequal(neighbor(1:end-1,2),(2:N)'));
mesh = struct('count',N,'points',p,'faces',b.faces,'originalCellId',order, ...
    'nodeIds',vertcat(nodes{order}),'boxes',boxes,'extent',extent, ...
    'neighbor',neighbor,'faceIds',faceIds,'sourceFile',cfg.meshFile);
mesh.dx = boxes(:,2)-boxes(:,1);
mesh.dy = boxes(:,4)-boxes(:,3);
assert(all(mesh.dx>tol) && all(mesh.dy>tol));
mesh.area = mesh.dx.*mesh.dy;
mesh.center = [(boxes(:,1)+boxes(:,2))/2,(boxes(:,3)+boxes(:,4))/2];
end

function vel = velocity_grid(order,temperature,wallModel)
if nargin < 3; wallModel = 'nonlinear'; end
[vx,wx] = hermite_rule(order(1));
[vy,wy] = hermite_rule(order(2));
[xx,yy] = meshgrid(vx,vy);
weights = wy*wx'/pi;
vel.vx = xx(:); vel.vy = yy(:); vel.weight = weights(:);
vel.speed2 = vel.vx.^2+vel.vy.^2;
vel.M0 = exp(-vel.speed2)/pi;
vel.basis = [ones(numel(xx),1),vel.vx,vel.vy,vel.speed2-1];
vel.moment = vel.weight.*vel.basis;
vel.gram = vel.basis'*vel.moment;
paired = flipud(reshape(1:numel(xx),order(2),order(1)));
vel.reflect = paired(:);
% Treat numerical roots at zero as exactly zero, so no arbitrary sweep is used.
vel.vx(abs(vel.vx)<1e-14) = 0;
vel.vy(abs(vel.vy)<1e-14) = 0;
vel.wallModel = wallModel;
vel.tauWall = temperature(:)' - 1;
vel.ratio = zeros(numel(xx),2);
for side = 1:2
    Tw = temperature(side);
    vel.ratio(:,side) = exp(vel.speed2*(1-1/Tw))/Tw;
end
assert(all(isfinite(vel.ratio),'all'),'Wall Maxwellian ratios overflow; change velocity order/temperature.');
assert(max(abs(vel.vx(vel.reflect)-vel.vx))<1e-12 && ...
    max(abs(vel.vy(vel.reflect)+vel.vy))<1e-12);
assert(max(abs(vel.weight(vel.reflect)-vel.weight))<1e-12);
vel.inLeft = find(vel.vx>0); vel.outLeft = find(vel.vx<0);
vel.inRight = vel.outLeft; vel.outRight = vel.inLeft;
vel.denominator = [sum(vel.weight(vel.inLeft).*vel.vx(vel.inLeft).*vel.ratio(vel.inLeft,1)), ...
    sum(vel.weight(vel.inRight).*(-vel.vx(vel.inRight)).*vel.ratio(vel.inRight,2))];
vel.linearDenominator = [sum(vel.weight(vel.inLeft).*vel.vx(vel.inLeft)), ...
    sum(vel.weight(vel.inRight).*(-vel.vx(vel.inRight)))];
vel.linearThermalShift = [ ...
    sum(vel.weight(vel.inLeft).*vel.vx(vel.inLeft).*(vel.speed2(vel.inLeft)-1))/vel.linearDenominator(1), ...
    sum(vel.weight(vel.inRight).*(-vel.vx(vel.inRight)).*(vel.speed2(vel.inRight)-1))/vel.linearDenominator(2)];
end

function [x,w] = hermite_rule(n)
% Integral exp(-x^2) p(x) dx, exact for polynomial degree <= 2*n-1.
a = sqrt((1:n-1)/2);
[U,D] = eig(diag(a,1)+diag(a,-1));
[x,idx] = sort(diag(D));
w = sqrt(pi)*U(1,idx)'.^2;
x = (x-flipud(x))/2;
w = (w+flipud(w))/2;
end

function op = operators(mesh,vel,delta)
N = mesh.count;
m = [2 1;1 2]/6;
d = [-1 -1;1 1]/2; % integral (test derivative)*(trial) on [0,1]
I = eye(4);
E = {I(:,[1 3]),I(:,[2 4]),I(:,[1 2]),I(:,[3 4])};
op.M = kron(spdiags(mesh.area,0,N,N),sparse(kron(m,m)));
op.massFactor = decomposition(op.M,'chol');
op.Vx = kron(spdiags(mesh.dy,0,N,N),sparse(-kron(m,d)));
op.Vy = kron(spdiags(mesh.dx,0,N,N),sparse(-kron(d,m)));
for side = 1:4
    if side<=2; length = mesh.dy; else; length = mesh.dx; end
    op.B{side} = kron(spdiags(length,0,N,N),sparse(E{side}*m*E{side}'));
end
op.Cleft = coupling(mesh,E{1}*m*E{2}',1);
op.Cright = coupling(mesh,E{2}*m*E{1}',2);
op.integral = repelem(mesh.area/4,4);
op.volume = sum(mesh.area);
op.transport = cell(numel(vel.vx),1);
op.factor = cell(numel(vel.vx),1);
for j = 1:numel(vel.vx)
    vx = vel.vx(j); vy = vel.vy(j);
    % Off-diagonal blocks are cross-cell matrices, never self-face matrices.
    T = delta*op.M+vx*op.Vx+vy*op.Vy + ...
        max(vx,0)*(op.B{2}-op.Cleft)+max(-vx,0)*(op.B{1}-op.Cright) + ...
        max(vy,0)*op.B{4}+max(-vy,0)*op.B{3};
    op.transport{j} = T;
    op.factor{j} = decomposition(T,'lu');
end
z = [-1;1]/sqrt(3);
op.edgeShape = [(1-z)/2,(1+z)/2];
op.edgeWeight = [0.5;0.5];
op.leftNodes = [1 3];
op.rightNodes = 4*(N-1)+[2 4];
end

function C = coupling(mesh,reference,side)
rows = []; cols = []; values = [];
[a,b,s] = find(reference);
for i = 1:mesh.count
    k = mesh.neighbor(i,side);
    if k>0
        rows = [rows;4*(i-1)+a]; %#ok<AGROW>
        cols = [cols;4*(k-1)+b]; %#ok<AGROW>
        values = [values;mesh.dy(i)*s]; %#ok<AGROW>
    end
end
C = sparse(rows,cols,values,4*mesh.count,4*mesh.count);
end

function [x,y] = nodal_coordinates(mesh)
b = mesh.boxes;
x = reshape([b(:,1),b(:,2),b(:,1),b(:,2)]',[],1);
y = reshape([b(:,3),b(:,3),b(:,4),b(:,4)]',[],1);
end

function H = field_values(value,x,y,vel)
if isa(value,'function_handle')
    H = value(x,y,vel.vx',vel.vy');
else
    H = value;
end
if isscalar(H); H = repmat(H,numel(x),numel(vel.vx)); end
assert(isequal(size(H),[numel(x),numel(vel.vx)]) && isreal(H) && all(isfinite(H),'all'), ...
    'initial/source must be a scalar or a function returning (4*Ncell)-by-Nvelocity.');
end

function [H,affine] = constrain_mass(H,op,vel,target,sourceFree)
rho = 1+(op.integral'*(H*vel.weight))/op.volume;
if sourceFree && ~strcmp(vel.wallModel,'linear')
    % Scaling the full f preserves the affine kinetic equation and diffuse
    % wall law when G=0, including unequal prescribed wall temperatures.
    assert(rho>0 && isfinite(rho),'BGKDG:NonpositiveMass','Nonpositive total full-distribution mass.');
    affine = [target/rho,target/rho-1];
else
    % For compatible nonzero sources or linearized diffuse walls, impose
    % the mass constraint as an additive shift (gauge).
    % Convergence still requires the ORIGINAL, unprojected kinetic residual.
    affine = [1,target-rho];
end
H = affine(1)*H+affine(2);
end

function wall = wall_data(H,mesh,vel,op)
shape = op.edgeShape;
left = shape*H(op.leftNodes,:);
right = shape*H(op.rightNodes,:);
leftOut = vel.outLeft; rightOut = vel.outRight;
if strcmp(vel.wallModel,'linear')
    nL = (left(:,leftOut)*(vel.weight(leftOut).*abs(vel.vx(leftOut))))/vel.linearDenominator(1) ...
        - vel.linearThermalShift(1)*vel.tauWall(1);
    nR = (right(:,rightOut)*(vel.weight(rightOut).*abs(vel.vx(rightOut))))/vel.linearDenominator(2) ...
        - vel.linearThermalShift(2)*vel.tauWall(2);
    rhoL = 1+nL;
    rhoR = 1+nR;
    gL = nL + vel.tauWall(1)*(vel.speed2-1)';
    gR = nR + vel.tauWall(2)*(vel.speed2-1)';
else
    rhoL = ((1+left(:,leftOut))*(vel.weight(leftOut).*abs(vel.vx(leftOut))))/vel.denominator(1);
    rhoR = ((1+right(:,rightOut))*(vel.weight(rightOut).*abs(vel.vx(rightOut))))/vel.denominator(2);
    gL = rhoL*vel.ratio(:,1)'-1;
    gR = rhoR*vel.ratio(:,2)'-1;
end
load = zeros(size(H));
incoming = vel.inLeft;
load(op.leftNodes,incoming) = mesh.dy(1)*(shape'*(op.edgeWeight.*gL(:,incoming))).*abs(vel.vx(incoming))';
incoming = vel.inRight;
load(op.rightNodes,incoming) = mesh.dy(end)*(shape'*(op.edgeWeight.*gR(:,incoming))).*abs(vel.vx(incoming))';
bottom = find(vel.vy>0); top = find(vel.vy<0);
load(:,bottom) = load(:,bottom)+(op.B{3}*H(:,vel.reflect(bottom))).*vel.vy(bottom)';
load(:,top) = load(:,top)+(op.B{4}*H(:,vel.reflect(top))).*(-vel.vy(top))';
wall = struct('load',load,'rhoLeft',rhoL,'rhoRight',rhoR, ...
    'gLeft',gL,'gRight',gR,'traceLeft',left,'traceRight',right, ...
    'reflectedH',H(:,vel.reflect));
end

function diagnostics = steady_diagnostics(H,Q,G,mesh,vel,op,cfg,usedWall)
% Recompute BOTH collision and wall inflows from the current H. This tests
% the steady equation, not merely the accuracy of a lagged local solve.
wall = wall_data(H,mesh,vel,op);
gain = (Q/vel.gram)*vel.basis';
R = -op.M*(cfg.delta*gain+G)-wall.load;
for j = 1:numel(vel.vx)
    R(:,j) = R(:,j)+op.transport{j}*H(:,j);
end
strong = op.massFactor\R;
diagnostics.kinetic = sqrt(max(0,sum(sum(R.*strong,1)'.*vel.weight)/op.volume));
diagnostics.massError = abs(1+(op.integral'*(H*vel.weight))/op.volume-cfg.meanDensity);
diagnostics.wallDensity = [wall.rhoLeft,wall.rhoRight];
% Compare the actual lagged inflows used by the last sweep to the physical
% wall map recomputed from CURRENT outgoing traces. DG does not require a
% strong match between the interior inflow polynomial and numerical trace.
[fluxLeft,fluxRight] = diffuse_mass_flux(wall,usedWall,vel);
diagnostics.diffuseMassFlux = [fluxLeft,fluxRight];
mass = vel.weight.*vel.vy;
bottom = H; top = H;
inB = find(vel.vy>0); inT = find(vel.vy<0);
bottom(:,inB) = usedWall.reflectedH(:,inB);
top(:,inT) = usedWall.reflectedH(:,inT);
bottomFlux = bottom*mass;
topFlux = top*mass;
bottomNodal = reshape(bottomFlux,4,mesh.count);
topNodal = reshape(topFlux,4,mesh.count);
diagnostics.specularMassFlux = [mean(bottomNodal(1:2,:),1)',mean(topNodal(3:4,:),1)'];
diagnostics.wallMassFlux = max(abs([fluxLeft(:);fluxRight(:);bottomFlux([1:4:end,2:4:end]);topFlux([3:4:end,4:4:end])]));
diffuseL = weighted_edge_error(wall.gLeft-usedWall.gLeft,vel.inLeft,vel,op);
diffuseR = weighted_edge_error(wall.gRight-usedWall.gRight,vel.inRight,vel,op);
change = wall.reflectedH-usedWall.reflectedH;
specB = boundary_change(change,inB,op.B{3},vel,op.volume);
specT = boundary_change(change,inT,op.B{4},vel,op.volume);
diagnostics.wall = max([diffuseL,diffuseR,specB,specT]);
diagnostics.wallComponents = [diffuseL,diffuseR,specB,specT];
diagnostics.collisionInvariantDefect = max(abs((H-gain)*vel.moment),[],'all');
diagnostics.totalMass = op.integral'*(1+H*vel.weight);
diagnostics.massSourceIntegral = op.integral'*(G*vel.weight);
end

function error = weighted_edge_error(change,ids,vel,op)
error = sqrt(sum((op.edgeWeight'*change(:,ids).^2)'.*vel.weight(ids)));
end

function error = boundary_change(change,ids,B,vel,volume)
% Norm of changed same-wall reflected outgoing traces, integrated along x.
error = sqrt(max(0,sum(sum(change(:,ids).*(B*change(:,ids)),1)'.*vel.weight(ids))/volume));
end

function [left,right] = diffuse_mass_flux(wall,usedWall,vel)
L = wall.traceLeft; R = wall.traceRight;
L(:,vel.inLeft) = usedWall.gLeft(:,vel.inLeft);
R(:,vel.inRight) = usedWall.gRight(:,vel.inRight);
mass = vel.weight.*vel.vx;
left = -(1+L)*mass;
right = (1+R)*mass;
end

function average = cell_average(nodal,N)
average = reshape(mean(reshape(nodal,4,N,[]),1),N,[]);
end

function plot_progress(window,mesh,Q,history)
figure(window); clf(window); tiledlayout(2,2);
q = cell_average(Q,mesh.count);
nexttile; plot(mesh.center(:,1),q(:,1),'LineWidth',1.2); xlabel('x'); ylabel('density perturbation');
nexttile; plot(mesh.center(:,1),q(:,4),'LineWidth',1.2); xlabel('x'); ylabel('temperature perturbation');
nexttile; plot(mesh.center(:,1),q(:,2:3),'LineWidth',1.2); xlabel('x'); ylabel('velocity moments'); legend('u_x','u_y');
nexttile; hold on;
for k = 1:3
    ids = isfinite(history(:,k));
    semilogy(history(ids,7),max(history(ids,k),1e-16),'LineWidth',1.1);
end
set(gca,'YScale','log'); xlabel('sweep'); ylabel('error / residual'); legend('moment update','kinetic','wall'); grid on;
drawnow;
end

function write_results(result)
folder = result.config.outputDir;
if ~exist(folder,'dir'); mkdir(folder); end
save(fullfile(folder,'bgk_dg.mat'),'result','-v7.3');
profile = [result.mesh.center,result.cellMoments,result.heatFlux];
writematrix(profile,fullfile(folder,'bgk_dg_profile.dat'),'Delimiter','tab');
writematrix(result.history,fullfile(folder,'bgk_dg_history.dat'),'Delimiter','tab');
fid = fopen(fullfile(folder,'bgk_dg_status.txt'),'w');
assert(fid>=0,'Unable to write solver status.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'converged=%d\niterations=%d\nwall_time=%.9g\n',result.converged,result.iterations,result.wallTime);
fprintf(fid,'kinetic_residual=%.12e\nwall_residual=%.12e\nwall_mass_flux=%.12e\nmass_error=%.12e\n', ...
    result.diagnostics.kinetic,result.diagnostics.wall,result.diagnostics.wallMassFlux,result.diagnostics.massError);
fprintf(fid,'profile columns: x y rho_pert ux uy T_pert heat_x heat_y\n');
fprintf(fid,'history columns: %s\n',strjoin(result.historyColumns,' '));
end
