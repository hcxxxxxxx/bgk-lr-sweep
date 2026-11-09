function check_sweep_structure
%CHECK_SWEEP_STRUCTURE Verify that T_j is block bidiagonal and that
%   decomposition(T,'lu')\rhs reproduces an explicit upwind cell sweep.
%   Rebuilds the same operators as monoBGKDG/operators on a uniform strip.
N = 40; dx = 0.25; dy = 0.25; delta = 0.3;
m = [2 1;1 2]/6;
d = [-1 -1;1 1]/2;
I = eye(4);
E = {I(:,[1 3]),I(:,[2 4]),I(:,[1 2]),I(:,[3 4])};
area = repmat(dx*dy,N,1);
M  = kron(spdiags(area,0,N,N),sparse(kron(m,m)));
Vx = kron(spdiags(repmat(dy,N,1),0,N,N),sparse(-kron(m,d)));
Vy = kron(spdiags(repmat(dx,N,1),0,N,N),sparse(-kron(d,m)));
B = cell(1,4);
for side = 1:4
    if side<=2; len = repmat(dy,N,1); else; len = repmat(dx,N,1); end
    B{side} = kron(spdiags(len,0,N,N),sparse(E{side}*m*E{side}'));
end
neighbor = zeros(N,4);
neighbor(2:end,1) = (1:N-1)';
neighbor(1:end-1,2) = (2:N)';
Cleft  = coupling(neighbor,E{1}*m*E{2}',1,dy,N);
Cright = coupling(neighbor,E{2}*m*E{1}',2,dy,N);

for trial = 1:4
    switch trial
        case 1; vx =  0.7; vy =  0.4;
        case 2; vx = -1.3; vy =  0.9;
        case 3; vx =  2.1; vy = -0.6;
        case 4; vx = -0.5; vy = -1.7;
    end
    T = delta*M+vx*Vx+vy*Vy + ...
        max(vx,0)*(B{2}-Cleft)+max(-vx,0)*(B{1}-Cright) + ...
        max(vy,0)*B{4}+max(-vy,0)*B{3};

    % (a) block structure: which 4x4 blocks are nonzero?
    pattern = zeros(N,N);
    for i = 1:N
        for k = 1:N
            blk = T(4*(i-1)+(1:4),4*(k-1)+(1:4));
            pattern(i,k) = nnz(blk)>0;
        end
    end
    [ri,ci] = find(pattern);
    offsets = unique(ci-ri)';
    if vx > 0
        expected = [-1 0];   % own block and upstream (left) block
    else
        expected = [0 1];
    end
    assert(isequal(sort(offsets),sort(expected)), ...
        'Unexpected block offsets %s for vx=%g.',mat2str(offsets),vx);

    % (b) fill-in of the cached factorisation
    [L,U,~,~] = lu(T);
    fillRatio = (nnz(L)+nnz(U))/nnz(T);

    % (c) explicit upwind sweep must match the cached backsolve
    rng(trial); rhs = randn(4*N,1);
    F = decomposition(T,'lu');
    hFactored = F\rhs;
    hSweep = explicit_sweep(T,rhs,N,vx);
    err = norm(hFactored-hSweep,inf)/norm(hFactored,inf);
    assert(err < 1e-10,'Sweep mismatch %.3e for vx=%g.',err,vx);

    % (d) diagonal block must be positive definite in its symmetric part
    D = full(T(1:4,1:4));
    lam = min(eig((D+D')/2));

    fprintf(['vx=%+5.2f vy=%+5.2f | block offsets %s | nnz(T)=%d ' ...
        'fill=%.2f | sweep err %.2e | min eig sym(D) %+.3e\n'], ...
        vx,vy,mat2str(sort(offsets)),nnz(T),fillRatio,err,lam);
end
fprintf('PASS: T is block bidiagonal and the cached solve equals the upwind sweep.\n');
end

function C = coupling(neighbor,reference,side,dy,N)
rows = []; cols = []; values = [];
[a,b,s] = find(reference);
for i = 1:N
    k = neighbor(i,side);
    if k>0
        rows = [rows;4*(i-1)+a]; %#ok<AGROW>
        cols = [cols;4*(k-1)+b]; %#ok<AGROW>
        values = [values;dy*s];  %#ok<AGROW>
    end
end
C = sparse(rows,cols,values,4*N,4*N);
end

function h = explicit_sweep(T,rhs,N,vx)
% Sequential cell-by-cell solve using the already-updated upstream block.
h = zeros(4*N,1);
if vx > 0; order = 1:N; up = -1; else; order = N:-1:1; up = +1; end
for i = order
    idx = 4*(i-1)+(1:4);
    load = rhs(idx);
    k = i+up;
    if k>=1 && k<=N
        kdx = 4*(k-1)+(1:4);
        load = load - T(idx,kdx)*h(kdx);
    end
    h(idx) = T(idx,idx)\load;
end
end
