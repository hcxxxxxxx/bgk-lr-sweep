function check_reference_stopping_1d2v
%CHECK_REFERENCE_STOPPING_1D2V Reference formula, source monitor and stopping.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
diary(fullfile(root,'debug','reference_stopping_validation.log'));
cleanup=onCleanup(@() diary('off'));
edges=[-.5;-.2;.5]; h=diff(edges);
Q=[.1 .03 .04 .2;.1 .03 .04 .2;-.2 0 .02 .05;-.2 0 .02 .05];
% Independently calculated cell integrals, including reference density/T.
expected=[.33 .36 .015;.56 .735 .014]; previous=[h,h,zeros(2,1)];
[metric,parts,current]=reference_macro_update_1d2v(Q,edges,previous);
expectedParts=sqrt(sum((expected-previous).^2,1)./sum(expected.^2,1));
assert(max(abs(current-expected),[],'all')<1e-14);
assert(max(abs(parts-expectedParts))<1e-14 && abs(metric-max(expectedParts))<1e-14);
flipped=Q; flipped(:,2:3)=-flipped(:,2:3);
assert(reference_macro_update_1d2v(flipped,edges,current)<1e-14);
assert(reference_macro_update_1d2v(zeros(4,4),edges)==0);
[~,~,moving]=reference_macro_update_1d2v(Q,edges);
assert(isinf(reference_macro_update_1d2v(zeros(4,4),edges,moving)));
fprintf('PASS reference formula on nonuniform cells, physical rho/T, speed magnitude and zero denominators.\n');

% A zero exact equilibrium can stop at once, without NaN speed ratios.
base=struct('flow','Couette','wallVelocity',[0 0],'cells',4, ...
    'normalQuadrature','hermite','velocityOrder',[8 8], ...
    'plot',false,'save',false,'reference',false);
z=mono1D2V(base);
assert(strcmp(z.config.stopCriterion,'reference') && z.config.tol==1e-5);
assert(z.converged && z.iterations==1 && z.convergence.referenceMacroUpdate==0);

% Fixed nonzero translating gas: first reference speed update is 1, second
% is zero, exactly as the Fortran OLD fields initialized to rest imply.
base.wallVelocity=[.007 .007]; base.initial=@(x,vx,vy) 0*x+0*vx+2*.007*vy;
uniform=monoGSIS1D2V(base);
assert(uniform.converged && uniform.iterations==2);
assert(abs(uniform.history(1,end)-1)<1e-14);

% GSIS must monitor returned source macros (as Calculate_Macro_Properties_ACC),
% not Qkin, and must not veto reference stopping using the equation norm.
opt=base; opt.flow='Poiseuille'; opt.wallVelocity=[0 0]; opt.initial=0;
opt.delta=10; opt.cells=8; opt.acceleration=1; opt.tol=.99; opt.maxIter=2;
g=monoGSIS1D2V(opt);
assert(g.converged && g.iterations==2 && g.diagnostics.coupled>opt.tol);
oldSource=g.synthetic.sourceMoments; oldSource(:,1)=oldSource(:,1)-mean(oldSource(:,1));
[~,~,sourceMonitor]=reference_macro_update_1d2v(oldSource,g.mesh.edges);
[~,~,kineticMonitor]=reference_macro_update_1d2v(g.nodalMoments,g.mesh.edges);
assert(max(abs(sourceMonitor-g.convergence.cellMonitor),[],'all')<1e-14);
assert(norm(sourceMonitor-kineticMonitor,'fro')>1e-6);
assert(isequal(g.synthetic.hdg.cache.stabilization,[1 10 10 10]));
fprintf('PASS source-macro monitoring and reference-only stopping, with unchanged stabilization.\n');
fprintf('All reference stopping checks passed.\n');
end
