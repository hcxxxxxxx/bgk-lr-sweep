function check_hdg_macro_1d2v
%CHECK_HDG_MACRO_1D2V Independent manufactured/monolithic/kinetic checks.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
diary(fullfile(root,'debug','hdg_macro_validation.log'));
cleanup=onCleanup(@() diary('off'));
mesh=struct('edges',[-.5;-.31;-.08;.17;.5]);
param=struct('stabilization',[.7 1.1 .4 .9]); delta=2.5;

% All four nonzero constants, including a nonzero normal velocity.
state=[.03 .004 .02 -.01]; bc=struct('trace',repmat(state,2,1)); % [p ux uy T]
primitiveState=state; primitiveState(1)=state(1)-state(4);
[q,a,h]=hdg_solve_macro_1d2v(mesh,delta,[0 0],0,0,0,bc,param);
assert(max(abs(q.nodal-primitiveState),[],'all')<1e-12);
assert(max(abs([a.Lx;a.Ly;a.E]))<1e-12);
check_residuals(h);
fprintf('PASS constant state on nonuniform cells, all seven fields.\n');
one=struct('edges',[-.5;.5]);
[single,~,hs]=hdg_solve_macro_1d2v(one,delta,[0 0],0,0,0,bc,param);
assert(max(abs(single.nodal-primitiveState),[],'all')<1e-12); check_residuals(hs);

% Exact P1 manufactured solution with BOTH accelerations and nonzero HOTs.
% ux=c; T=t0+t1*x; Pi_xx=p0+p1*x; uy=u0+u1*x;
% Pi_xy=b0-2*delta*ay*x gives (sigma_xy/2)'=ay.
% p is chosen so [p-Lx/delta]'=2*ax.
x=h.cache.nodes; ax=.013; ay=.007;
Pi_xx=.01+.03*x; Pi_xy=.02-2*delta*ay*x; Theta=.007+0*x;
exact=manufactured(x,delta,ax);
bc=@(ctx) struct('trace',manufactured(ctx.boundaryX,delta,ax));
[q,a,h]=hdg_solve_macro_1d2v(mesh,delta,[ax ay],Pi_xx,Pi_xy,Theta,bc,param);
assert(max(abs(q.pressureNodal-exact),[],'all')<1e-12);
assert(max(abs(q.rho-(exact(:,1)-exact(:,4))))<1e-12);
assert(max(abs(a.Lx-Pi_xx))<1e-12);
assert(max(abs(a.Ly-(.04+Pi_xy)))<1e-12);
assert(max(abs(a.E-.017))<1e-12);
for i=1:numel(mesh.edges)-1
    ids=2*i-1:2*i;
    assert(abs(diff(a.Pxy(ids))/diff(mesh.edges(i:i+1))-ay)<1e-12);
end
check_residuals(h);
fprintf('PASS all-field affine manufactured solution, normal forcing and Pxy''=ay.\n');

% Identical data, cached vs uncached operator.
cached=param; cached.cache=h.cache;
[again,~,hc]=hdg_solve_macro_1d2v(mesh,delta,[ax ay],Pi_xx,Pi_xy,Theta,bc,cached);
assert(max(abs(again.nodal-q.nodal),[],'all')<1e-13);
assert(max(abs(hc.Uhat-h.Uhat),[],'all')<1e-13);
% A cache also carries non-default stabilization when only .cache is passed.
[again,~,~]=hdg_solve_macro_1d2v(mesh,delta,[ax ay], ...
    @(z) .01+.03*z,@(z) .02-2*delta*ay*z,.007,bc,struct('cache',h.cache));
assert(max(abs(again.nodal-q.nodal),[],'all')<1e-13);

% Monolithic system contains every cell field AND every face trace.
% Uses the original weak operator blocks, without the Schur complement.
bcMixed=struct('trace',[.01 .002 -.01 .003;.02 .002 NaN .006], ...
    'flux',[NaN NaN NaN NaN;NaN NaN .004 NaN]);
N=numel(mesh.edges)-1; forcing=[sin((1:2*N)'),cos((1:2*N)')]/100;
[~,~,h]=hdg_solve_macro_1d2v(mesh,delta,forcing,Pi_xx,Pi_xy,Theta,bcMixed,param);
[fullLocal,fullTrace]=monolithic(h);
assert(max(abs(fullLocal-h.localFields),[],'all')<1e-11);
assert(max(abs(fullTrace-h.Uhat),[],'all')<1e-11);
check_residuals(h);
fprintf('PASS static condensation vs monolithic mixed-boundary solve: field %.3e, trace %.3e.\n', ...
    max(abs(fullLocal-h.localFields),[],'all'),max(abs(fullTrace-h.Uhat),[],'all'));

% Couette with zero HOT has the exact affine solution and constant stress.
mesh=struct('edges',linspace(-.5,.5,17)');
bc=struct('trace',[0 0 -.01 0;0 0 .01 0]);
[q,a,h]=hdg_solve_macro_1d2v(mesh,10,[0 0],0,0,0,bc);
assert(isequal(h.cache.stabilization,[1 10 10 10]));
assert(max(abs(q.uy-.02*h.cache.nodes))<1e-12);
assert(max(abs(a.Pxy+.001))<1e-12); check_residuals(h);
fprintf('PASS Couette analytic profile and constant Pxy.\n');
for strength=[.01 1000]
    [qs,as,hs]=hdg_solve_macro_1d2v(mesh,strength,[0 0],0,0,0,bc);
    assert(isequal(hs.cache.stabilization,[1 strength strength strength]));
    assert(max(abs(qs.uy-.02*hs.cache.nodes))<1e-10);
    assert(max(abs(as.Pxy+.01/strength))<1e-10);
    check_residuals(hs);
end
fprintf('PASS macroscopic Couette operator at delta=0.01 and 1000.\n');

% Poiseuille: numerical flux balance at machine precision and spatial
% convergence to u=ay*delta*(1/4-x^2), with zero HOT/no-slip BC.
bc=struct('trace',zeros(2,4)); errors=zeros(3,1);
for k=1:3
    N=16*2^(k-1); mesh=struct('edges',linspace(-.5,.5,N+1)');
    [q,~,h]=hdg_solve_macro_1d2v(mesh,3,[0 .01],0,0,0,bc);
    exact=.03*(.25-h.cache.nodes.^2);
    errors(k)=sqrt(sum((q.uy-exact).^2)/(2*N));
    assert(max(abs(h.outwardFlux(3,:)+h.outwardFlux(7,:)-.02/N))<1e-12);
    check_residuals(h);
end
fprintf('Poiseuille RMS errors N=16/32/64: %.6e %.6e %.6e\n',errors);
assert(all(errors(1:2)./errors(2:3)>3.5));
fprintf('PASS Poiseuille conservative force balance and second-order spatial convergence.\n');

% Get Gaussian nodes from the existing kinetic implementation.
k=mono1D2V(struct('flow','Couette','cells',4,'normalQuadrature','hermite', ...
    'velocityOrder',[8 8],'plot',false,'save',false,'reference',false, ...
    'stopCriterion','equation','tol',1e-10));
v=k.velocity; w=v.weight;
assert(abs(sum(w)-1)<1e-13 && abs(sum(w.*v.vx.^2)-.5)<1e-13);
assert(abs(sum(w.*v.vx.^2.*v.vy.^2)-.25)<1e-13);
assert(isequal(size(k.kineticTrace),[5 numel(w)]));
equilibrium=primitiveState*v.basis';
hot=compute_high_order_moments_1d2v(repmat(equilibrium,8,1),v);
assert(max(abs(hot.R),[],'all')<1e-13);
mesh=k.mesh;
face=compute_high_order_moments_1d2v(repmat(equilibrium,5,1),v);
d=project_hot_derivatives_1d(mesh,hot.R,face.R);
assert(max(abs(d),[],'all')<1e-12);
H3=permute(reshape(repmat(equilibrium,8,1),2,4,[]),[2 1 3]);
hot3=compute_high_order_moments_1d2v(H3,v);
assert(isequal(hot3.R,hot.R));
fprintf('PASS velocity normalization, 2D Gaussian moments, equilibrium HOT and H layouts.\n');

% A DG jump with zero volume Rxy but NONZERO upwind face Rxy.
% Taking a moment first and averaging would incorrectly give zero.
mesh=struct('edges',[-.5;0;.5]);
H=[1;1;3;3]*(v.vx.*v.vy)';
trace=zeros(3,numel(w)); trace(1,:)=H(1,:); trace(3,:)=H(4,:);
trace(2,v.vx>0)=H(2,v.vx>0); trace(2,v.vx<0)=H(3,v.vx<0);
hot=compute_high_order_moments_1d2v(H,v);
face=compute_high_order_moments_1d2v(trace,v);
weak=project_hot_derivatives_1d(mesh,hot.R,face.R);
strong=project_hot_derivatives_1d(mesh,hot.R,[],'strong');
assert(max(abs(hot.Rxy))<1e-13 && abs(face.Rxy(2))>.01);
assert(max(abs(weak(:,2)))>.01 && max(abs(strong(:,2)))<1e-12);
for i=1:2
    M=.5*[2 1;1 2]/6;
    assert(max(abs(sum(M*weak(2*i-1:2*i,:),1)-(face.R(i+1,:)-face.R(i,:))))<1e-12);
end
fprintf('PASS velocity-by-velocity upwind trace, DG jumps and weak-derivative integral identity.\n');

% Actual current 32x32 Couette/Poiseuille -> HOT, kinetic BC, HDG subsolve.
cases=mono1D2V(struct('plot',false,'save',false,'reference',false, ...
    'stopCriterion','equation','tol',1e-10));
for name={'couette','poiseuille'}
    k=cases.(name{1});
    [q,~,h,hot]=gsis_macro_from_kinetic_1d2v(k);
    check_residuals(h);
    assert(max(abs(hot.boundary.trace(:,2)))<1e-12);
    assert(max(abs(h.Uhat([1 end],:)-hot.boundary.trace),[],'all')<1e-12);
    assert(strcmp(hot.derivativeMethod,'interior'));
    expectedBoundary=[k.nodalMoments([1 end],1)+k.nodalMoments([1 end],4), ...
        zeros(2,1),k.nodalMoments([1 end],3:4)];
    assert(max(abs(hot.boundary.trace-expectedBoundary),[],'all')<1e-12);
    change=norm(q.nodal(:,3)-k.nodalMoments(:,3))/norm(k.nodalMoments(:,3));
    fprintf('PASS kinetic-to-HDG %s: local %.3e, trace %.3e, kinetic uy change %.6e\n', ...
        name{1},h.localResidual,h.traceResidual,change);
end
fprintf('All 1D2V HDG macroscopic checks passed.\n');
end

function values=manufactured(x,delta,ax)
T=.005+.01*x;
values=[.02+(2*ax+.03/delta)*x,.003+0*x,.01+.04*x,T];
end

function check_residuals(h)
assert(max([h.traceResidual,h.localResidual,h.fluxContinuityResidual,h.conservationResidual])<1e-11);
end

function [U,hat]=monolithic(h)
c=h.cache; N=c.count; nf=14*N; nt=4*(N+1);
% This independent reference system is tiny; dense assembly avoids relying
% on the production sparse trace assembly or its elimination order.
A=zeros(nf+nt,nf+nt); b=zeros(nf+nt,1);
for i=1:N
    fields=14*(i-1)+(1:14); traces=nf+c.traceIds(:,i);
    A(fields,fields)=c.A(:,:,i); A(fields,traces)=-c.B(:,:,i);
    A(traces,fields)=A(traces,fields)+c.C(:,:,i);
    A(traces,traces)=A(traces,traces)+c.D(:,:,i);
    b(fields)=h.localRhs(:,i);
end
% Prescribed natural fluxes are nonzero only on boundary rows.
for side=1:2
    ids=(side-1)*4*N+(1:4);
    values=h.bc.flux(side,:); active=isfinite(values);
    b(nf+ids(active))=values(active);
end
for j=1:numel(h.bc.known)
    row=nf+h.bc.known(j); A(row,:)=0; A(row,row)=1; b(row)=h.bc.knownValues(j);
end
solution=A\b;
assert(max(abs(A*solution-b))<1e-11);
U=reshape(solution(1:nf),14,N); hat=reshape(solution(nf+1:end),4,[])';
end
