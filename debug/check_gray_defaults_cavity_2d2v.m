function check_gray_defaults_cavity_2d2v
% Exercise reference defaults, native-style stopping, and unchanged H.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
old=warning('off','Cavity:NotConverged');cleanup=onCleanup(@()warning(old));
for delta=[3 1000]
    cfg=struct('cells',[3 2],'meshStretch',1.3,'velocityOrder',[8 8], ...
        'quadrature','hermite','delta',delta,'maxIter',1,'plot',false,'save',false);
    cis=monoCavity2D2V(cfg);g=monoGSISCavity2D2V(cfg);
    assert(isequal(cis.H,g.H));
    assert(g.config.tol==1e-5 && ~g.synthetic.hdg.internalResidualCorrection);
    assert(g.synthetic.relaxation==1 && isequal(g.synthetic.hdg.stabilization,ones(1,4)));
    assert(~any(g.synthetic.hdg.loadCorrection,'all'));
    assert(g.synthetic.hdg.hotUsesTransportWallInflow);
    assert(strcmp(g.convergenceCriterion,'referenceCellIntegratedMacroUpdate'));
    assert(max(abs(g.synthetic.sourceMoments-g.synthetic.macro.nodal),[],'all')<1e-12);
    [expected,components]=native_loop(g.synthetic.sourceMoments,zeros(size(g.nodalMoments)),g.mesh,g.config.tol);
    assert(abs(g.history(end,2)-expected)<1e-14);
    assert(max(abs(g.history(end,10:13)-components))<1e-14);
    assert(g.converged==(expected<g.config.tol));
end
% The two wall data paths must stay distinct: changing the stored last-sweep
% inflow changes HOT by exactly its face moment, but not the affine wall load.
c=hdg_precompute_q1_2d2v(cis.mesh,cis.velocity,cis.config.delta);
[~,base]=gsis_macro_from_kinetic_2d2v(cis,struct('cache',c));
changed=cis;cellId=1;side=1;nodes=[1 3];increment=[.003;-.002];
changed.transportWallDensity(nodes,side)=changed.transportWallDensity(nodes,side)+increment;
[~,altered]=gsis_macro_from_kinetic_2d2v(changed,struct('cache',c));
v=cis.velocity;x=v.vx;y=v.vy;s=v.speed2;cn=-x;incoming=cn<0;
rx=[x.*(x.^2-y.^2-1),y.*(2*x.^2-1),x.^2.*(s-2)-(s-1),x.*y.*(s-2)];
normalMoment=-sum(v.weight(incoming).*rx(incoming,:),1);
expected=c.E{side}*c.faceM(:,:,side,cellId)*(increment*normalMoment);
actual=altered.hotWeakLoad(:,:,cellId)-base.hotWeakLoad(:,:,cellId);
assert(max(abs(actual-expected),[],'all')<1e-13 && norm(actual,'fro')>1e-8);
assert(isequal(altered.wallAffineTarget,base.wallAffineTarget));
assert(isequal(altered.wallKineticFlux,base.wallKineticFlux));
% Independently reproduce the archive's cell-integral loop and both sides
% of its small-denominator branch, including a nearly zero component.
mesh=g.mesh;Q=zeros(4*mesh.count,4);previous=Q;
Q(:,1)=.1*mesh.nodes(:,1);Q(:,2)=100*mesh.nodes(:,2);
Q(:,3)=1e-12;Q(:,4)=.2*mesh.nodes(:,2);
previous(:,2)=Q(:,2)-1e-3;previous(:,3)=-1e-12;
for tol=[1e-5 1e-8]
    [actual,component]=reference_macro_update_2d2v(Q,previous,mesh,tol);
    [expected,components]=native_loop(Q,previous,mesh,tol);
    assert(abs(actual-expected)<1e-14 && max(abs(component-components))<1e-14);
    assert(component(3)<1e-10,'Near-zero component must use reference-norm floor.');
end
cfg=struct('cells',4,'velocityOrder',[8 8],'quadrature','hermite','delta',10, ...
    'tol',1e-5,'maxIter',300,'plot',false,'save',false);
g=monoGSISCavity2D2V(cfg);
assert(g.converged && g.history(end,2)<cfg.tol);
assert(max(g.history(end,10:13))==g.history(end,2));
assert(strcmp(g.synthetic.hdg.formulation,'gray-bgk-2v-v2'));
fprintf('PASS Gray-aligned 2D defaults at delta=3/1000, no H lift or extra load, unit stabilization, full return, cell-integral stopping and small-component floor.\n');
fprintf('Reference-default smoke: delta=10, steps=%d, macro update %.6e, original relative %.6e, coupled relative %.6e.\n', ...
    g.iterations,g.history(end,2),g.diagnostics.relativeKinetic,g.diagnostics.relativeCoupled);
end

function [residual,component]=native_loop(Q,previous,mesh,tol)
ss1=zeros(1,4);ss2=zeros(1,4);
for cellId=1:mesh.count
    nodes=4*(cellId-1)+(1:4);
    current=mesh.area(cellId)*(sum(Q(nodes,:),1)/4+[1 0 0 1]);
    old=mesh.area(cellId)*(sum(previous(nodes,:),1)/4+[1 0 0 1]);
    ss1=ss1+(current-old).^2;ss2=ss2+current.^2;
end
reference=max(max(ss2),1);component=zeros(1,4);
for k=1:4
    if ss2(k)>tol*tol*reference;component(k)=sqrt(ss1(k)/ss2(k));
    else;component(k)=sqrt(ss1(k)/reference);end
end
residual=max(component);
end
