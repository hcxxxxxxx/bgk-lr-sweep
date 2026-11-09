function bench_sweep
%BENCH_SWEEP Cost of the cached backsolve versus alternatives.
%   Compares, for one velocity system T_j h = r:
%     (1) decomposition(T,'lu') built once, then backsolve   <- monoBGKDG
%     (2) T\r  (fresh factorisation every call)
%     (3) explicit cell-by-cell upwind sweep in MATLAB
%     (4) lu(T) factors reused with two triangular solves
N = 400; dx = 0.25; dy = 0.25; delta = 0.01;
vx = 0.7; vy = 0.4;
T = build_transport(N,dx,dy,delta,vx,vy);
rng(0); rhs = randn(4*N,1);

F = decomposition(T,'lu');
ref = F\rhs;

% (1) cached decomposition backsolve
t1 = timeit(@() F\rhs);

% (2) fresh factorisation each solve
t2 = timeit(@() T\rhs);

% (3) explicit sweep written in MATLAB
t3 = timeit(@() explicit_sweep(T,rhs,N,vx));

% (4) lu factors + two triangular solves
[L,U,P,Qp] = lu(T);
t4 = timeit(@() Qp*(U\(L\(P*rhs))));

% (5) cost of building the decomposition itself
t5 = timeit(@() decomposition(T,'lu'));

err3 = norm(explicit_sweep(T,rhs,N,vx)-ref,inf)/norm(ref,inf);
err4 = norm(Qp*(U\(L\(P*rhs)))-ref,inf)/norm(ref,inf);

fprintf('N = %d cells, 4N = %d unknowns, nnz(T) = %d\n\n',N,4*N,nnz(T));
fprintf('%-38s %12s %10s %10s\n','method','time [us]','ratio','rel err');
fprintf('%-38s %12.2f %10.2f %10.1e\n','(1) cached decomposition backsolve',t1*1e6,1.0,0);
fprintf('%-38s %12.2f %10.2f %10.1e\n','(2) T\\rhs, refactor every call',t2*1e6,t2/t1,0);
fprintf('%-38s %12.2f %10.2f %10.1e\n','(3) explicit MATLAB cell sweep',t3*1e6,t3/t1,err3);
fprintf('%-38s %12.2f %10.2f %10.1e\n','(4) reused lu + triangular solves',t4*1e6,t4/t1,err4);
fprintf('%-38s %12.2f %10.2f %10s\n','(5) building decomposition (once)',t5*1e6,t5/t1,'-');

nvel = 256; iters = 6989;
fprintf('\nProjected cost for %d velocities x %d sweeps:\n',nvel,iters);
fprintf('  (1) cached backsolve      : %8.1f s   (+ %.2f s setup)\n', ...
    t1*nvel*iters, t5*nvel);
fprintf('  (2) refactor every sweep  : %8.1f s\n',t2*nvel*iters);
fprintf('  (3) explicit MATLAB sweep : %8.1f s\n',t3*nvel*iters);
end

function T = build_transport(N,dx,dy,delta,vx,vy)
m = [2 1;1 2]/6;
d = [-1 -1;1 1]/2;
I = eye(4);
E = {I(:,[1 3]),I(:,[2 4]),I(:,[1 2]),I(:,[3 4])};
M  = kron(spdiags(repmat(dx*dy,N,1),0,N,N),sparse(kron(m,m)));
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
T = delta*M+vx*Vx+vy*Vy + ...
    max(vx,0)*(B{2}-Cleft)+max(-vx,0)*(B{1}-Cright) + ...
    max(vy,0)*B{4}+max(-vy,0)*B{3};
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
