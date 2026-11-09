function bench_sweep2
%BENCH_SWEEP2 Loop-based timing, closer to how monoBGKDG actually calls it.
N = 400; dx = 0.25; dy = 0.25; delta = 0.01;
vx = 0.7; vy = 0.4;
T = build_transport(N,dx,dy,delta,vx,vy);
rng(0); R = randn(4*N,64);          % 64 right-hand sides
reps = 64;

F = decomposition(T,'lu');
[L,U,P,Qp] = lu(T);

% warm up
tmp = F\R(:,1); tmp = U\(L\(P*R(:,1))); %#ok<NASGU>

t = zeros(1,5);

tic; for k = 1:reps; x = F\R(:,k); end; t(1) = toc/reps;
tic; for k = 1:reps; x = T\R(:,k); end; t(2) = toc/reps;
tic; for k = 1:reps; x = Qp*(U\(L\(P*R(:,k)))); end; t(3) = toc/reps;
tic; for k = 1:reps; x = explicit_sweep(T,R(:,k),N,vx); end; t(4) = toc/reps;
tic; for k = 1:reps; Fk = decomposition(T,'lu'); end; t(5) = toc/reps; %#ok<NASGU>

% block solve: all 64 rhs at once (what a vectorised code could do)
tic; X = F\R; tB = toc;

fprintf('N=%d, 4N=%d unknowns, nnz(T)=%d, %d right-hand sides\n\n',N,4*N,nnz(T),reps);
names = {'cached decomposition  F\\r', ...
         'T\\r (refactor each time)', ...
         'reused lu: Q*(U\\(L\\(P*r)))', ...
         'explicit MATLAB cell sweep', ...
         'decomposition(T,''lu'') build'};
fprintf('%-34s %12s %8s\n','method','time [us]','ratio');
for i = 1:5
    fprintf('%-34s %12.2f %8.2f\n',names{i},t(i)*1e6,t(i)/t(1));
end
fprintf('%-34s %12.2f %8.2f\n','F\\R all 64 columns at once',tB/reps*1e6,tB/reps/t(1));

nvel = 256; iters = 6989;
fprintf('\nProjected for %d velocities x %d sweeps (solve only):\n',nvel,iters);
for i = [1 2 3 4]
    fprintf('  %-32s %8.1f s\n',names{i},t(i)*nvel*iters);
end
fprintf('  %-32s %8.1f s\n','F\\R batched over velocities',tB/reps*nvel*iters);
fprintf('\nActual monoBGKDG run: 147.7 s total for the same workload.\n');
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
