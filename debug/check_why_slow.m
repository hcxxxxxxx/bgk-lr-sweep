function check_why_slow
%CHECK_WHY_SLOW Diagnose why monoBGKDG takes ~3582 steps to reach 1e-5 at delta=0.01
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);

r = monoBGKDG(struct('plot',false,'save',false,'maxIter',1,'velocityOrder',[16 16],'delta',0.01));
mesh = r.mesh; vel = r.velocity; cfg = r.config;
N = mesh.count; V = numel(vel.vx); nDof = 4*N;

% Build operators
m = [2 1;1 2]/6; d = [-1 -1;1 1]/2; I = eye(4);
E = {I(:,[1 3]),I(:,[2 4]),I(:,[1 2]),I(:,[3 4])};
M = kron(spdiags(mesh.area,0,N,N),sparse(kron(m,m)));
massFactor = decomposition(M,'chol');
Vx = kron(spdiags(mesh.dy,0,N,N),sparse(-kron(m,d)));
Vy = kron(spdiags(mesh.dx,0,N,N),sparse(-kron(d,m)));
B = cell(4,1);
for side = 1:4
    if side<=2; len = mesh.dy; else; len = mesh.dx; end
    B{side} = kron(spdiags(len,0,N,N),sparse(E{side}*m*E{side}'));
end
Cleft  = build_coupling(mesh,E{1}*m*E{2}',1);
Cright = build_coupling(mesh,E{2}*m*E{1}',2);
integral = repelem(mesh.area/4,4);
volume = sum(mesh.area);
z = [-1;1]/sqrt(3);
shape = [(1-z)/2,(1+z)/2];
edgeWeight = [0.5;0.5];
leftNodes = [1 3]; rightNodes = 4*(N-1)+[2 4];

% Notice: for this 1D channel problem, the exact solution is symmetric in y
% and even in vy: H(:, reflect(j)) == H(:, j)!
% In monoBGKDG, top/bottom specular reflection is LAGGED:
%   LHS has +max(vy,0)*B{4} + max(-vy,0)*B{3}
%   RHS has +max(vy,0)*B{3}*H_old(:,jr) + max(-vy,0)*B{4}*H_old(:,jr)
% Whereas if specular reflection is UNLAGGED (coupled (j,jr) or even-vy):
%   For an even-in-vy field H(:,jr) = H(:,j), the top/bottom inflow term moves
%   directly to the LHS: -max(vy,0)*B{3}*H(:,j) - max(-vy,0)*B{4}*H(:,j)!

T_orig = cell(V,1); F_orig = cell(V,1);
T_unlag_spec = cell(V,1); F_unlag_spec = cell(V,1);
for j = 1:V
    vx = vel.vx(j); vy = vel.vy(j);
    T = cfg.delta*M + vx*Vx + vy*Vy + ...
        max(vx,0)*(B{2}-Cleft) + max(-vx,0)*(B{1}-Cright) + ...
        max(vy,0)*B{4} + max(-vy,0)*B{3};
    T_orig{j} = T;
    F_orig{j} = decomposition(T,'lu');
    % Unlagged specular for vy-even solution (H(:,reflect(j)) = H(:,j)):
    T2 = T - max(vy,0)*B{3} - max(-vy,0)*B{4};
    T_unlag_spec{j} = T2;
    F_unlag_spec{j} = decomposition(T2,'lu');
end

fprintf('\n====================================================================\n');
fprintf('  Comparing Iteration Counts to tol = 1e-5 (delta = 0.01, Lx = 100)\n');
fprintf('====================================================================\n');

for mode = 1:3
    H = zeros(nDof,V);
    Q = H*vel.moment;
    for iter = 1:500
        prevQ = Q;
        gain = (Q/vel.gram)*vel.basis';
        if mode == 1
            % Mode 1: Original monoBGKDG (lagged specular + lagged diffuse)
            left = shape*H(leftNodes,:); right = shape*H(rightNodes,:);
            nL = (left(:,vel.outLeft)*(vel.weight(vel.outLeft).*abs(vel.vx(vel.outLeft))))/vel.linearDenominator(1) ...
                - vel.linearThermalShift(1)*vel.tauWall(1);
            nR = (right(:,vel.outRight)*(vel.weight(vel.outRight).*abs(vel.vx(vel.outRight))))/vel.linearDenominator(2) ...
                - vel.linearThermalShift(2)*vel.tauWall(2);
            gL = nL + vel.tauWall(1)*(vel.speed2-1)';
            gR = nR + vel.tauWall(2)*(vel.speed2-1)';
            rhs = cfg.delta*(M*gain);
            in = vel.inLeft;
            rhs(leftNodes,in) = rhs(leftNodes,in) + mesh.dy(1)*(shape'*(edgeWeight.*gL(:,in))).*abs(vel.vx(in))';
            in = vel.inRight;
            rhs(rightNodes,in) = rhs(rightNodes,in) + mesh.dy(end)*(shape'*(edgeWeight.*gR(:,in))).*abs(vel.vx(in))';
            bot = find(vel.vy>0); top = find(vel.vy<0);
            rhs(:,bot) = rhs(:,bot) + (B{3}*H(:,vel.reflect(bot))).*vel.vy(bot)';
            rhs(:,top) = rhs(:,top) + (B{4}*H(:,vel.reflect(top))).*(-vel.vy(top))';
            for j = 1:V
                H(:,j) = F_orig{j}\rhs(:,j);
            end
        elseif mode == 2
            % Mode 2: Unlagged specular + Lagged diffuse
            left = shape*H(leftNodes,:); right = shape*H(rightNodes,:);
            nL = (left(:,vel.outLeft)*(vel.weight(vel.outLeft).*abs(vel.vx(vel.outLeft))))/vel.linearDenominator(1) ...
                - vel.linearThermalShift(1)*vel.tauWall(1);
            nR = (right(:,vel.outRight)*(vel.weight(vel.outRight).*abs(vel.vx(vel.outRight))))/vel.linearDenominator(2) ...
                - vel.linearThermalShift(2)*vel.tauWall(2);
            gL = nL + vel.tauWall(1)*(vel.speed2-1)';
            gR = nR + vel.tauWall(2)*(vel.speed2-1)';
            rhs = cfg.delta*(M*gain);
            in = vel.inLeft;
            rhs(leftNodes,in) = rhs(leftNodes,in) + mesh.dy(1)*(shape'*(edgeWeight.*gL(:,in))).*abs(vel.vx(in))';
            in = vel.inRight;
            rhs(rightNodes,in) = rhs(rightNodes,in) + mesh.dy(end)*(shape'*(edgeWeight.*gR(:,in))).*abs(vel.vx(in))';
            for j = 1:V
                H(:,j) = F_unlag_spec{j}\rhs(:,j);
            end
        elseif mode == 3
            % Mode 3: Unlagged specular + Unlagged (Gauss-Seidel) diffuse
            rhs = cfg.delta*(M*gain);
            right = shape*H(rightNodes,:);
            nR = (right(:,vel.outRight)*(vel.weight(vel.outRight).*abs(vel.vx(vel.outRight))))/vel.linearDenominator(2) ...
                - vel.linearThermalShift(2)*vel.tauWall(2);
            gR = nR + vel.tauWall(2)*(vel.speed2-1)';
            in = vel.inRight;
            rhs(rightNodes,in) = rhs(rightNodes,in) + mesh.dy(end)*(shape'*(edgeWeight.*gR(:,in))).*abs(vel.vx(in))';
            for j = in'
                H(:,j) = F_unlag_spec{j}\rhs(:,j);
            end
            left = shape*H(leftNodes,:);
            nL = (left(:,vel.outLeft)*(vel.weight(vel.outLeft).*abs(vel.vx(vel.outLeft))))/vel.linearDenominator(1) ...
                - vel.linearThermalShift(1)*vel.tauWall(1);
            gL = nL + vel.tauWall(1)*(vel.speed2-1)';
            in = vel.inLeft;
            rhs(leftNodes,in) = rhs(leftNodes,in) + mesh.dy(1)*(shape'*(edgeWeight.*gL(:,in))).*abs(vel.vx(in))';
            for j = in'
                H(:,j) = F_unlag_spec{j}\rhs(:,j);
            end
        end
        rho = 1 + (integral'*(H*vel.weight))/volume;
        H = H + (1 - rho);
        Q = H*vel.moment;
        update = max(sqrt(sum((Q-prevQ).^2,1)./max(sum(prevQ.^2,1),size(Q,1))));
        % Compute true steady kinetic residual
        left = shape*H(leftNodes,:); right = shape*H(rightNodes,:);
        nL = (left(:,vel.outLeft)*(vel.weight(vel.outLeft).*abs(vel.vx(vel.outLeft))))/vel.linearDenominator(1) ...
            - vel.linearThermalShift(1)*vel.tauWall(1);
        nR = (right(:,vel.outRight)*(vel.weight(vel.outRight).*abs(vel.vx(vel.outRight))))/vel.linearDenominator(2) ...
            - vel.linearThermalShift(2)*vel.tauWall(2);
        gL = nL + vel.tauWall(1)*(vel.speed2-1)';
        gR = nR + vel.tauWall(2)*(vel.speed2-1)';
        R = -cfg.delta*(M*((Q/vel.gram)*vel.basis'));
        in = vel.inLeft;
        R(leftNodes,in) = R(leftNodes,in) - mesh.dy(1)*(shape'*(edgeWeight.*gL(:,in))).*abs(vel.vx(in))';
        in = vel.inRight;
        R(rightNodes,in) = R(rightNodes,in) - mesh.dy(end)*(shape'*(edgeWeight.*gR(:,in))).*abs(vel.vx(in))';
        for j = 1:V
            R(:,j) = R(:,j) + T_unlag_spec{j}*H(:,j);
        end
        strong = massFactor\R;
        kinetic = sqrt(max(0,sum(sum(R.*strong,1)'.*vel.weight)/volume));
        if mode > 1 && (iter <= 15 || kinetic < 1e-5)
            fprintf('  Mode %d, step %2d: momentUpdate = %.3e, kineticResidual = %.3e\n', mode, iter, update, kinetic);
        end
        if kinetic < 1e-5 && update < 1e-5
            fprintf('--> Mode %d CONVERGED to < 1e-5 in %d steps!\n\n', mode, iter);
            break;
        end
    end
end
end

function C = build_coupling(mesh,reference,side)
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
