function summary=run_quintic_cavity_delta1000(n)
% Quintic wall-clustered Q1 grid; preserve Gray source return and S=I.
% Continue beyond native stopping so weak corner velocities can settle.
if nargin<1;n=64;end
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
folder=fullfile(root,'results','cavity_delta1000_quintic_20261112');
if ~exist(folder,'dir');mkdir(folder);end
diary(fullfile(folder,sprintf('run_N%d.log',n)));cleanup=onCleanup(@()diary('off'));
cfg=struct('cells',n,'meshMapping','quintic','delta',1000, ...
    'velocityOrder',[32 32],'tol',1e-5,'tolCoupled',1e-9,'minIter',60,'maxIter',400, ...
    'reportEvery',10,'plot',false,'outputDir',folder);
r=monoGSISCavity2D2V(cfg);
assert(r.converged && r.iterations>=60 && all(isfinite(r.H),'all'));
assert(r.diagnostics.relativeCoupled<cfg.tolCoupled);
assert(isequal(r.synthetic.hdg.stabilization,ones(1,4)));
assert(~r.synthetic.hdg.internalResidualCorrection && r.synthetic.relaxation==1);
assert(all(diff(r.mesh.xEdges)>0) && max(abs(r.mesh.xEdges+fliplr(r.mesh.xEdges)-1))<1e-14);
caseFolder=fullfile(folder,sprintf('delta_1000_N%dx%d_quintic_GSIS',n,n));
report=plotCavityCorners2D2V(r,caseFolder);
c=report.ellipticVelocityRoots;
first=find(r.history(:,2)<r.config.tol,1);
U=reshape(r.cellMoments(:,2),n,n);V=reshape(r.cellMoments(:,3),n,n);
symmetry=max([max(abs(U-flipud(U)),[],'all'),max(abs(V+flipud(V)),[],'all')]);
summary=table(n,min(diff(r.mesh.xEdges)),max(diff(r.mesh.xEdges)), ...
    first,r.iterations,r.diagnostics.macroUpdate,r.diagnostics.relativeKinetic, ...
    r.diagnostics.relativeCoupled,sum(c(:,1)<.2 & c(:,2)<.2), ...
    sum(c(:,1)>.8 & c(:,2)<.2),symmetry,r.wallTime, ...
    'VariableNames',{'cellsPerDirection','minimumCellWidth','maximumCellWidth', ...
    'firstNativeStop','iterations','macroUpdate','originalRelativeResidual', ...
    'coupledRelativeResidual','leftCornerCandidates','rightCornerCandidates', ...
    'velocitySymmetryError','seconds'});
writetable(summary,fullfile(folder,sprintf('summary_N%d.csv',n)));disp(summary);
end
