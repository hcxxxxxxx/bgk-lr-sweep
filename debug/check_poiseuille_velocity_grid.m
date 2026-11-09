% Isolate velocity quadrature effects in the published-Qp comparison.
% Fix Nx=256, gap=1, A=1, stationary diffuse walls, tol=1e-10.
% Reuse the saved Nx256, Vx512, Vy8, cutoff8 run as the baseline.
% Compare Vx256/512/1024, cutoff8/10, and Vy8/16 ONE setting at a time.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
folder = fullfile(root,'results','channel_1d2v','literature_poiseuille');
baseline = readtable(fullfile(folder,'comparison.csv'));
diary(fullfile(folder,'velocity_check.log'));
diaryCleanup = onCleanup(@() diary('off'));
cfg = struct('flow','Poiseuille','gap',1,'acceleration',1, ...
    'wallVelocity',[0 0],'initial',0,'tol',1e-10,'stopCriterion','equation','maxIter',10000, ...
    'cells',256,'normalQuadrature','halfLegendre', ...
    'checkEvery',10,'reportEvery',1000,'reference',false,'plot',false,'save',false);
% Rows: Nvx, Nvy, cutoff. Middle row is reused from the previous comparison.
variants = [256 8 8; 512 8 8; 1024 8 8; 1024 8 10; 512 16 8];
rows = zeros(height(baseline)*size(variants,1),12);
n = 0;
for k = 1:height(baseline)
    cfg.delta = baseline.delta(k);
    for j = 1:size(variants,1)
        cfg.velocityOrder = variants(j,1:2); cfg.normalCutoff = variants(j,3);
        if j==2
            qp = baseline.fineQp(k); iterations = baseline.fineIterations(k);
            residual = baseline.fineResidual(k);
        else
            r = mono1D2V(cfg);
            assert(r.converged,'Velocity check did not converge.');
            assert(max(abs(r.nodalMoments(:,[1 2 4])),[],'all')<1e-9);
            assert(max(abs(r.velocity.vx))<cfg.normalCutoff);
            assert(abs(sum(r.velocity.weight)-1)<1e-12);
            qp = r.flowRate/(cfg.acceleration*cfg.gap^2);
            iterations = r.iterations; residual = r.diagnostics.kinetic;
        end
        n = n+1;
        rows(n,:) = [cfg.delta,cfg.cells,cfg.velocityOrder,cfg.normalCutoff, ...
            cfg.tol,qp,baseline.paperBGK(k), ...
            (qp/baseline.fineQp(k)-1),(qp/baseline.paperBGK(k)-1), ...
            iterations,residual];
        fprintf('VELOCITY CHECK: delta=%g Nv=%dx%d cutoff=%g Qp=%.14g change=%+.5e paper_diff=%+.5e\n', ...
            cfg.delta,cfg.velocityOrder,cfg.normalCutoff,qp,rows(n,9:10));
    end
end
velocityComparison = array2table(rows,'VariableNames', ...
    {'delta','cells','Nvx','Nvy','normalCutoff','tolerance','Qp','paperBGK', ...
    'relativeChangeFromBaseline','relativeDifferenceToPaper','iterations','kineticResidual'});
writetable(velocityComparison,fullfile(folder,'velocity_check.csv'));

% Report independent changes as relative fractions (not percentages).
summaryRows = zeros(height(baseline),5);
for k = 1:height(baseline)
    ids = (k-1)*size(variants,1)+(1:size(variants,1));
    q = velocityComparison.Qp(ids);
    summaryRows(k,:) = [baseline.delta(k),abs(q(2)/q(1)-1), ...
        abs(q(3)/q(2)-1),abs(q(4)/q(3)-1),abs(q(5)/q(2)-1)];
end
velocitySummary = array2table(summaryRows,'VariableNames', ...
    {'delta','Nvx256to512','Nvx512to1024','cutoff8to10_atNvx1024','Nvy8to16_atNvx512'});
writetable(velocitySummary,fullfile(folder,'velocity_check_summary.csv'));
save(fullfile(folder,'velocity_check.mat'),'velocityComparison','velocitySummary','cfg','variants');
disp(velocitySummary);
fprintf('All summary changes are relative fractions, not percentages.\n');
fprintf('Velocity quadrature check saved in %s\n',folder);
