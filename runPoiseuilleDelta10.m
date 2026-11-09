function r = runPoiseuilleDelta10
%RUNPOISEUILLEDELTA10 Compare delta=1 and 10 at the SAME acceleration 0.01.
%   r = runPoiseuilleDelta10;
% Results: results/channel_1d2v/delta10/poiseuille. The saved delta=1 run
% and delta=10 Couette run are retained. A finer DG run checks resolution.
root = fileparts(mfilename('fullpath'));
folder = fullfile(root,'results','channel_1d2v','delta10','poiseuille');
if ~exist(folder,'dir'); mkdir(folder); end
fid = fopen(fullfile(folder,'run.log'),'w'); assert(fid>=0); fclose(fid);
diary(fullfile(folder,'run.log')); cleanup = onCleanup(@() diary('off'));
cfg = struct('flow','Poiseuille','delta',1,'acceleration',0.01, ...
    'wallVelocity',[0 0],'cells',128,'velocityOrder',[128 8], ...
    'referenceCells',2048,'referenceOrder',320,'plot',false,'save',false, ...
    'reportEvery',200,'outputDir',fullfile(root,'results','channel_1d2v','delta10'));
baseline = mono1D2V(cfg);
cfg.delta = 10; cfg.plot = true; cfg.save = true;
r = mono1D2V(cfg);

% Check space and normal-speed refinement without overwriting the main run.
fineCfg = r.config; fineCfg.cells = 256; fineCfg.velocityOrder = [256 8];
fineCfg.plot = false; fineCfg.save = false; fineCfg.reference = false;
fine = mono1D2V(fineCfg);
fineReference = mean(reshape(r.reference.velocity,r.reference.cells/fineCfg.cells,fineCfg.cells),1)';
fineError = norm(fine.cellMoments(:,3)-fineReference)/norm(fineReference);
gridChange = norm(mean(reshape(fine.cellMoments(:,3),2,[]),1)'-r.cellMoments(:,3))/norm(r.cellMoments(:,3));
coarseReference = channel1D2VReference(r.config,1024,320);
refChange = norm(mean(reshape(r.reference.velocity,2,[]),1)'-coarseReference.velocity)/norm(coarseReference.velocity);
assert(baseline.converged && r.converged && fine.converged);
assert(r.reference.relativeL2<1e-4 && fineError<1e-4 && gridChange<1e-4 && refChange<1e-4);
assert(max(abs(r.nodalMoments(:,[1 2 4])),[],'all')<1e-12);
assert(max(r.diagnostics.balance)<1e-8);
assert(max(abs(r.slipVelocity))>0 && all(r.cellMoments(:,3)>0));
assert(max(abs(r.cellMoments(:,3)-flipud(r.cellMoments(:,3))))<1e-12);

fprintf('\nComparison: acceleration=0.01, stationary diffuse walls, gap=1\n');
fprintf(' delta    max cell-mean uy      flow rate        wall uy           NS flow rate\n');
metrics = zeros(2,7);
for k = 1:2
    if k==1; s = baseline; else; s = r; end
    nsFlow = s.config.acceleration*s.config.delta*s.config.gap^3/6;
    peak = max(s.cellMoments(:,3));
    metrics(k,:) = [s.config.delta,peak,s.flowRate,s.gasWallVelocity(1),nsFlow, ...
        s.reference.relativeL2,s.iterations];
    fprintf('%6g    %.10g       %.10g     %.10g      %.10g\n',metrics(k,1:5));
    fprintf('  wall/peak fraction %.6g, flow excess over no-slip NS %.6g\n', ...
        s.gasWallVelocity(1)/peak,(s.flowRate-nsFlow)/nsFlow);
end
fprintf('delta=10 DG relative L2 vs BGK reference: N128/Vx128 %.5e, N256/Vx256 %.5e\n', ...
    r.reference.relativeL2,fineError);
fprintf('DG refinement change %.5e; reference N1024->N2048 change %.5e\n',gridChange,refChange);
fprintf('max |ux,rho'',T''| %.5e; max invariant balance defect %.5e\n', ...
    max(abs(r.nodalMoments(:,[1 2 4])),[],'all'),max(r.diagnostics.balance));

window = figure('Name','Poiseuille delta 1 vs 10','Color','w','Position',[100 100 1150 780]);
if isprop(window,'Theme'); window.Theme='light'; end
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
baselineX = [-0.5;baseline.mesh.center;0.5];
baselineU = [baseline.gasWallVelocity(1);baseline.cellMoments(:,3);baseline.gasWallVelocity(2)];
currentX = [-0.5;r.mesh.center;0.5];
currentU = [r.gasWallVelocity(1);r.cellMoments(:,3);r.gasWallVelocity(2)];
plot(baselineX,baselineU,'LineWidth',1.4);
plot(currentX,currentU,'LineWidth',1.4);
plot(r.mesh.center,r.reference.cellVelocity,'k--','LineWidth',1.1);
xx = linspace(-0.5,0.5,501)';
plot(xx,0.01*10*(0.25-xx.^2),':','Color',[.4 .4 .4],'LineWidth',1.2);
xlabel('x'); ylabel('u_y'); grid on;
legend('delta=1 DG','delta=10 DG','delta=10 BGK reference','delta=10 NS no slip','Location','best');
title('Same acceleration: a=0.01');

nexttile; hold on;
plot(baselineX,baselineU/max(baseline.cellMoments(:,3)),'LineWidth',1.3);
plot(currentX,currentU/max(r.cellMoments(:,3)),'LineWidth',1.3);
plot(xx,1-4*xx.^2,'k:','LineWidth',1.2);
xlabel('x'); ylabel('u_y / max cell-mean(u_y)'); grid on;
legend('delta=1','delta=10','NS no slip shape','Location','best');
title('Profile shape approaches the no-slip parabola');

nexttile; hold on;
plot(baseline.mesh.center,baseline.shearStress,'LineWidth',1.3);
plot(r.mesh.center,r.shearStress,'--','LineWidth',1.3);
plot(xx,0.01*xx,'k:','LineWidth',1.2);
xlabel('x'); ylabel('tau_{xy}'); grid on;
legend('delta=1','delta=10','a x','Location','best');
title('Both satisfy d(tau_{xy})/dx = a');

nexttile; hold on;
ids = isfinite(baseline.history(:,3));
semilogy(baseline.history(ids,1),baseline.history(ids,3),'LineWidth',1.3);
ids = isfinite(r.history(:,3));
semilogy(r.history(ids,1),r.history(ids,3),'LineWidth',1.3);
set(gca,'YScale','log'); xlabel('source iteration'); ylabel('kinetic residual'); grid on;
legend('delta=1','delta=10','Location','best'); title('Ordinary source iteration slows as delta grows');
sgtitle('Poiseuille: delta=1 and 10, same force and wall conditions'); drawnow;
exportgraphics(window,fullfile(folder,'delta_comparison.png'),'Resolution',150);
writematrix(metrics,fullfile(folder,'delta_comparison.dat'),'Delimiter','tab');
fid = fopen(fullfile(folder,'delta_comparison.txt'),'w'); assert(fid>=0);
fileCleanup = onCleanup(@() fclose(fid));
fprintf(fid,'Same acceleration a=0.01, wallVelocity=[0 0], gap=1.\n');
fprintf(fid,'columns: delta max_cell_mean_uy flow_rate left_wall_uy NS_no_slip_flow_rate BGK_reference_relative_L2 iterations\n');
fprintf(fid,'DG refinement relative change=%.12e\nfine_DG_reference_relative_L2=%.12e\nreference_spatial_refinement_relative_change=%.12e\n', ...
    gridChange,fineError,refChange);
fprintf(fid,'max_abs_ux_rhoPert_TPert=%.12e\nmax_invariant_balance_defect=%.12e\n', ...
    max(abs(r.nodalMoments(:,[1 2 4])),[],'all'),max(r.diagnostics.balance));
fprintf('Poiseuille delta=10 checks passed; output: %s\n',folder);
end
