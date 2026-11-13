function report=analyze_cavity_delta100_cis
% Inspect only the user-requested delta=100 pure CIS, tol=1e-5 case.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
folder=fullfile(root,'results','cavity_delta100_cis_tol1e5','delta_100_N32x32_stretch2');
data=load(fullfile(folder,'cavity.mat'),'r');r=data.r;
assert(strcmp(r.config.iterationMethod,'CIS') && r.config.delta==100 && r.config.tol==1e-5);
assert(r.converged && r.diagnostics.relativeKinetic<1e-5);
report=plotCavityCorners2D2V(r,folder);
roots=report.ellipticVelocityRoots;
bottom=roots(:,2)<.2 & (roots(:,1)<.2 | roots(:,1)>.8);
fprintf('Bottom-corner elliptic velocity roots: %d\n',sum(bottom));
disp(roots(bottom,:));
minWidth=min([diff(r.mesh.xEdges),diff(r.mesh.yEdges)]);
summary=table(r.config.delta,r.config.tol,r.iterations,r.diagnostics.macroUpdate, ...
    r.diagnostics.kinetic,r.diagnostics.relativeKinetic,r.diagnostics.wallMassFlux, ...
    minWidth,report.maximumEdgeFluxProjectionError,sum(bottom), ...
    'VariableNames',{'delta','tolerance','cisIterations','macroRelativeUpdate', ...
    'absoluteKineticResidual','relativeKineticResidual','lastSweepWallMassFlux', ...
    'minimumCellWidth','streamfunctionProjectionError','bottomCornerCandidates'});
writetable(summary,fullfile(folder,'corner_summary.csv'));disp(summary);
end
