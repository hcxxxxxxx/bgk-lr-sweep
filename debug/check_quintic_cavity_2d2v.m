function check_quintic_cavity_2d2v
% Verify unchanged defaults and kinetic solutions on the new mapped grid.
check_cavity_2d2v;
check_gray_defaults_cavity_2d2v;
c=struct('cells',8,'meshMapping','quintic','lidSpeed',0, ...
    'velocityOrder',[4 4],'quadrature','hermite','plot',false,'save',false, ...
    'maxIter',3,'minIter',3);
z=monoGSISCavity2D2V(c);
assert(z.converged && z.iterations==3 && max(abs(z.H),[],'all')==0);
c.lidSpeed=.01;c.cells=4;c.tol=1e-9;c.maxIter=1000;c.minIter=0;
c.iterationMethod='CIS';a=monoCavity2D2V(c);
c.iterationMethod='GMRES';b=monoCavity2D2V(c);
assert(a.converged && b.converged && max(abs(a.H-b.H),[],'all')<1e-9);
fprintf('PASS quintic grid equilibrium, minimum iterations, CIS vs GMRES.\n');
c=struct('cells',4,'meshMapping','quintic','delta',10, ...
    'velocityOrder',[8 8],'quadrature','hermite','plot',false,'save',false, ...
    'maxIter',400,'tol',1e-5);
native=monoGSISCavity2D2V(c);
c.tolCoupled=1e-10;strict=monoGSISCavity2D2V(c);
assert(strict.converged && strict.diagnostics.relativeCoupled<c.tolCoupled);
assert(strict.iterations>native.iterations);
assert(strcmp(strict.convergenceCriterion,'referenceMacroUpdateAndCoupledEquation'));
fprintf('PASS optional coupled-residual stopping extends iteration without changing defaults.\n');
end
