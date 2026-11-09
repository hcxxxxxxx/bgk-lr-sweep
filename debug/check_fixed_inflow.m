function check_fixed_inflow
%CHECK_FIXED_INFLOW Verify 1-step convergence under prescribed inflow BCs.
%   Compares three collisionless (delta -> 0) scenarios on mesh/nsw400.cas:
%     Case 1: Prescribed Dirichlet inflow on all 4 boundaries (L, R, B, T).
%     Case 2: Prescribed L/R inflow + LAGGED top/bottom specular reflection.
%     Case 3: Prescribed L/R inflow + UNLAGGED (implicit) top/bottom specular reflection.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);

% Grab mesh and velocity grid from monoBGKDG (1 dummy step, no plot/save)
r = monoBGKDG(struct('plot',false,'save',false,'maxIter',1,'velocityOrder',[8 8]));
mesh = r.mesh;
vel  = r.velocity;
N    = mesh.count;
V    = numel(vel.vx);
nDof = 4 * N;

% Assemble the exact same Q1-DG operators as monoBGKDG with delta = 0 (collisionless)
delta = 0;
op = build_operators(mesh, vel, delta);

% Define fixed left and right inflow distributions gL(v), gR(v)
% Example: linearized thermal wall perturbation g_w = tau_w * (v^2 - 1)
tauL = -0.05;
tauR =  0.05;
gL = repmat(tauL * (vel.speed2 - 1)', 2, 1); % 2 Gauss points on left edge x V
gR = repmat(tauR * (vel.speed2 - 1)', 2, 1); % 2 Gauss points on right edge x V

% Build fixed left/right inflow RHS load
rhsLR = zeros(nDof, V);
inL = vel.inLeft;
rhsLR(op.leftNodes, inL) = mesh.dy(1) * (op.edgeShape' * (op.edgeWeight .* gL(:, inL))) .* abs(vel.vx(inL))';
inR = vel.inRight;
rhsLR(op.rightNodes, inR) = mesh.dy(end) * (op.edgeShape' * (op.edgeWeight .* gR(:, inR))) .* abs(vel.vx(inR))';

fprintf('=========================================================================\n');
fprintf('  Collisionless Transport (delta = 0), N = %d cells, V = %d velocities\n', N, V);
fprintf('  Fixed Left/Right Inflow: tauL = %.2f, tauR = %.2f, Initial H^(0) = 0\n', tauL, tauR);
fprintf('=========================================================================\n\n');

%% ------------------------------------------------------------------------
%  CASE 1: All 4 boundaries have prescribed fixed inflow (Dirichlet L/R/B/T)
%  For a channel with y-independent inflow:
%    h*(x,y,v) = gL(v) for vx > 0, and gR(v) for vx < 0.
%  Suppose top/bottom inflow is also prescribed to match this boundary data.
%% ------------------------------------------------------------------------
H_exact_dirichlet = zeros(nDof, V);
H_exact_dirichlet(:, inL) = repmat(gL(1, inL), nDof, 1);
H_exact_dirichlet(:, inR) = repmat(gR(1, inR), nDof, 1);

rhs_all_fixed = rhsLR;
bottom = find(vel.vy > 0);
top    = find(vel.vy < 0);
rhs_all_fixed(:, bottom) = rhs_all_fixed(:, bottom) + (op.B{3} * H_exact_dirichlet(:, bottom)) .* vel.vy(bottom)';
rhs_all_fixed(:, top)    = rhs_all_fixed(:, top)    + (op.B{4} * H_exact_dirichlet(:, top))    .* (-vel.vy(top))';

H_case1 = zeros(nDof, V);
for j = 1:V
    H_case1(:, j) = op.factor{j} \ rhs_all_fixed(:, j);
end

res_case1 = zeros(nDof, V);
for j = 1:V
    res_case1(:, j) = op.transport{j} * H_case1(:, j) - rhs_all_fixed(:, j);
end
err_case1 = max(abs(H_case1 - H_exact_dirichlet), [], 'all');
norm_res1 = max(abs(res_case1), [], 'all');

fprintf('[Case 1] Fixed Inflow on ALL 4 Boundaries (L, R, B, T):\n');
fprintf('  After 1 sweep -> max |T*H - rhs|         = %.3e\n', norm_res1);
fprintf('  After 1 sweep -> max |H^(1) - H_exact|   = %.3e  (1-step exact!)\n\n', err_case1);

%% ------------------------------------------------------------------------
%  CASE 2: Fixed L/R Inflow + LAGGED Top/Bottom Specular Reflection
%  Start from H^(0) = 0. In sweep 1, top/bottom specular load uses H^(0) = 0.
%% ------------------------------------------------------------------------
H_lag = zeros(nDof, V); % H^(0) = 0
fprintf('[Case 2] Fixed L/R Inflow + LAGGED Top/Bottom Specular Reflection:\n');
for step = [1, 2, 5, 10, 50, 200]
    while true
        % Perform 1 lagged sweep from current H_lag
        rhs_lag = rhsLR;
        rhs_lag(:, bottom) = rhs_lag(:, bottom) + (op.B{3} * H_lag(:, vel.reflect(bottom))) .* vel.vy(bottom)';
        rhs_lag(:, top)    = rhs_lag(:, top)    + (op.B{4} * H_lag(:, vel.reflect(top)))    .* (-vel.vy(top))';
        H_next = zeros(nDof, V);
        for j = 1:V
            H_next(:, j) = op.factor{j} \ rhs_lag(:, j);
        end
        H_lag = H_next;
        if ~exist('curStep','var'); curStep = 1; else; curStep = curStep + 1; end
        if curStep == step; break; end
    end
    steady_res = compute_specular_steady_residual(H_lag, rhsLR, op, vel);
    err_exact  = max(abs(H_lag - H_exact_dirichlet), [], 'all');
    fprintf('  After %3d sweep(s) -> steady residual = %.3e,  max |H - H_exact| = %.3e\n', ...
        step, steady_res, err_exact);
end
fprintf('  --> Because Lx/Ly = 100/0.25 = 400, lagging top/bottom reflection wipes out\n');
fprintf('      the inflow with old H^(n) values every time a particle bounces on y=0 or y=0.25!\n\n');

%% ------------------------------------------------------------------------
%  CASE 3: Fixed L/R Inflow + UNLAGGED (Implicit) Top/Bottom Specular Reflection
%  Couple (j, reflect(j)) implicitly on top/bottom walls and solve in 1 sweep.
%% ------------------------------------------------------------------------
H_case3 = zeros(nDof, V);
visited = false(V, 1);
for j = 1:V
    if visited(j); continue; end
    jr = vel.reflect(j);
    visited([j, jr]) = true;
    if j == jr
        H_case3(:, j) = op.factor{j} \ rhsLR(:, j);
    else
        vy = vel.vy(j);
        if vy > 0
            % j has vy > 0 (bottom inflow B{3}), jr has vy < 0 (top inflow B{4})
            S12 = -vy * op.B{3};
            S21 = -vy * op.B{4};
        else
            % j has vy < 0 (top inflow B{4}), jr has vy > 0 (bottom inflow B{3})
            S12 =  vy * op.B{4};
            S21 =  vy * op.B{3};
        end
        T_pair = [op.transport{j}, S12; S21, op.transport{jr}];
        sol = T_pair \ [rhsLR(:, j); rhsLR(:, jr)];
        H_case3(:, j)  = sol(1:nDof);
        H_case3(:, jr) = sol(nDof+1:end);
    end
end

steady_res3 = compute_specular_steady_residual(H_case3, rhsLR, op, vel);
err_case3   = max(abs(H_case3 - H_exact_dirichlet), [], 'all');

fprintf('[Case 3] Fixed L/R Inflow + UNLAGGED (Implicit) Top/Bottom Specular Reflection:\n');
fprintf('  After 1 sweep -> steady residual         = %.3e\n', steady_res3);
fprintf('  After 1 sweep -> max |H^(1) - H_exact|   = %.3e  (1-step exact!)\n\n', err_case3);

%% ------------------------------------------------------------------------
%  BONUS CASE 4: Diffuse Walls (Linear Model, Unlagged) with delta = 0
%  Even without prescribing n_L and n_R a priori, because for delta = 0
%  left-going particles come purely from the right wall and right-going
%  particles come purely from the left wall, solving the 2x2 wall density
%  balance (or 2 sweeps) gives the exact collisionless solution immediately!
%% ------------------------------------------------------------------------
% In 1st pass: set nR^(0) = 0 -> sweep vx < 0 -> get exact outgoing flux at Left
% -> compute nL -> sweep vx > 0 -> get exact outgoing flux at Right -> compute nR!
% Moreover, for symmetric geometry with zero net mass (nL + nR = 0), nL = 0, nR = 0
% is already the exact steady state!
fprintf('=========================================================================\n');
end

function res_norm = compute_specular_steady_residual(H, rhsLR, op, vel)
% Evaluate the true unlagged steady residual:
%   R_j = T_j * H_j - rhsLR_j - specular_load(H_{reflect(j)})
R = -rhsLR;
bottom = find(vel.vy > 0);
top    = find(vel.vy < 0);
R(:, bottom) = R(:, bottom) - (op.B{3} * H(:, vel.reflect(bottom))) .* vel.vy(bottom)';
R(:, top)    = R(:, top)    - (op.B{4} * H(:, vel.reflect(top)))    .* (-vel.vy(top))';
for j = 1:numel(vel.vx)
    R(:, j) = R(:, j) + op.transport{j} * H(:, j);
end
strong = op.massFactor \ R;
res_norm = sqrt(max(0, sum(sum(R .* strong, 1)' .* vel.weight) / op.volume));
end

function op = build_operators(mesh, vel, delta)
N = mesh.count;
m = [2 1; 1 2] / 6;
d = [-1 -1; 1 1] / 2;
I = eye(4);
E = {I(:, [1 3]), I(:, [2 4]), I(:, [1 2]), I(:, [3 4])};
op.M = kron(spdiags(mesh.area, 0, N, N), sparse(kron(m, m)));
op.massFactor = decomposition(op.M, 'chol');
op.Vx = kron(spdiags(mesh.dy, 0, N, N), sparse(-kron(m, d)));
op.Vy = kron(spdiags(mesh.dx, 0, N, N), sparse(-kron(d, m)));
for side = 1:4
    if side <= 2; len = mesh.dy; else; len = mesh.dx; end
    op.B{side} = kron(spdiags(len, 0, N, N), sparse(E{side} * m * E{side}'));
end
op.Cleft  = build_coupling(mesh, E{1} * m * E{2}', 1);
op.Cright = build_coupling(mesh, E{2} * m * E{1}', 2);
op.volume = sum(mesh.area);
op.transport = cell(numel(vel.vx), 1);
op.factor    = cell(numel(vel.vx), 1);
for j = 1:numel(vel.vx)
    vx = vel.vx(j); vy = vel.vy(j);
    T = delta * op.M + vx * op.Vx + vy * op.Vy + ...
        max(vx, 0) * (op.B{2} - op.Cleft) + max(-vx, 0) * (op.B{1} - op.Cright) + ...
        max(vy, 0) * op.B{4} + max(-vy, 0) * op.B{3};
    op.transport{j} = T;
    op.factor{j}    = decomposition(T, 'lu');
end
z = [-1; 1] / sqrt(3);
op.edgeShape  = [(1 - z) / 2, (1 + z) / 2];
op.edgeWeight = [0.5; 0.5];
op.leftNodes  = [1 3];
op.rightNodes = 4 * (N - 1) + [2 4];
end

function C = build_coupling(mesh, reference, side)
rows = []; cols = []; values = [];
[a, b, s] = find(reference);
for i = 1:mesh.count
    k = mesh.neighbor(i, side);
    if k > 0
        rows   = [rows; 4 * (i - 1) + a]; %#ok<AGROW>
        cols   = [cols; 4 * (k - 1) + b]; %#ok<AGROW>
        values = [values; mesh.dy(i) * s]; %#ok<AGROW>
    end
end
C = sparse(rows, cols, values, 4 * mesh.count, 4 * mesh.count);
end
