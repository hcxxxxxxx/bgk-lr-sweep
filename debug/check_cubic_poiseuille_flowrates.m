function summary = check_cubic_poiseuille_flowrates
% All six published parallel-plate Poiseuille Qp values at fixed 32x32.
% Main comparison uses the agreed relative-macro stop 1e-5. A separate
% CIS-only 1e-9 check quantifies stopping sensitivity without mixing tables.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
folder=fullfile(root,'results','cubic_poiseuille_flowrates32');
if ~exist(folder,'dir');mkdir(folder);end
diary(fullfile(folder,'run.log'));cleanup=onCleanup(@()diary('off'));
delta=[.1 .5 1 2 5 10];
paperBGK=[2.0327 1.6022 1.5392 1.5957 1.9922 2.7673];
cfg=struct('flow','Poiseuille','gap',1,'cells',64,'acceleration',.01, ...
    'wallVelocity',[0 0],'initial',0,'velocityOrder',[32 32], ...
    'stopCriterion','reference','tol',1e-5,'maxIter',30000, ...
    'plot',false,'save',false,'reference',false,'reportEvery',1000);
rows=cell(0,13);strictRows=cell(0,7);diagnosticRows=cell(0,11);
runs=struct('delta',{},'grid',{},'cis',{},'gsis',{},'strictCIS',{});
for i=1:numel(delta)
    cfg.delta=delta(i);
    for j=1:2
        if j==1
            grid='original';cfg.normalQuadrature='halfLegendre';cfg.normalCutoff=8;
            cfg.tangentialQuadrature='hermite';cfg.tangentialCutoff=4;
        else
            grid='referenceCubic';cfg.normalQuadrature='referenceCubic';cfg.normalCutoff=4;
            cfg.tangentialQuadrature='referenceCubic';cfg.tangentialCutoff=4;
        end
        c=mono1D2V(cfg);g=monoGSIS1D2V(cfg);
        for method={'cis','gsis'}
            if strcmp(method{1},'cis');r=c;else;r=g;end
            check_run(r);
            diagnosticRows(end+1,:)={delta(i),grid,method{1},r.config.tol,r.iterations, ...
                r.convergence.referenceMacroUpdate,r.diagnostics.kinetic,r.diagnostics.coupled, ...
                r.diagnostics.collisionInvariantDefect,r.diagnostics.wallMassFlux,r.diagnostics.massError}; %#ok<AGROW>
        end
        qp=[c.flowRate g.flowRate]/(cfg.acceleration*cfg.gap^2);
        difference=100*(qp/paperBGK(i)-1);
        rows(end+1,:)={delta(i),grid,c.flowRate,g.flowRate,qp(1),qp(2),paperBGK(i), ...
            difference(1),difference(2),c.iterations,g.iterations, ...
            100*(qp(2)/qp(1)-1),c.iterations/g.iterations}; %#ok<AGROW>
        tight=[];
        if j==2
            tightCfg=cfg;tightCfg.tol=1e-9;tight=mono1D2V(tightCfg);check_run(tight);
            tightQp=tight.flowRate/(cfg.acceleration*cfg.gap^2);
            strictRows(end+1,:)={delta(i),tightQp,paperBGK(i), ...
                100*(tightQp/paperBGK(i)-1),tight.iterations, ...
                100*(tightQp-qp(1))/tightQp,tight.convergence.referenceMacroUpdate}; %#ok<AGROW>
            diagnosticRows(end+1,:)={delta(i),grid,'cis',tight.config.tol,tight.iterations, ...
                tight.convergence.referenceMacroUpdate,tight.diagnostics.kinetic,tight.diagnostics.coupled, ...
                tight.diagnostics.collisionInvariantDefect,tight.diagnostics.wallMassFlux,tight.diagnostics.massError}; %#ok<AGROW>
        end
        runs(end+1)=struct('delta',delta(i),'grid',grid,'cis',c,'gsis',g,'strictCIS',tight); %#ok<AGROW>
        fprintf('FLOW delta=%g %s: CIS Qp=%.12g (%+.8g%%), GSIS Qp=%.12g (%+.8g%%), steps %d/%d\n', ...
            delta(i),grid,qp(1),difference(1),qp(2),difference(2),c.iterations,g.iterations);
    end
end
summary=cell2table(rows,'VariableNames',{'delta','grid','cisQ','gsisQ','cisQp','gsisQp', ...
    'paperBGKQp','cisVsPaperPercent','gsisVsPaperPercent','cisIterations','gsisIterations', ...
    'signedGSISvsCISFlowPercent','iterationRatio'});
strictCIS=cell2table(strictRows,'VariableNames',{'delta','cisQp','paperBGKQp', ...
    'cisVsPaperPercent','cisIterations','signedStoppingEffectPercent','referenceMacroUpdate'});
diagnostics=cell2table(diagnosticRows,'VariableNames',{'delta','grid','method','tolerance', ...
    'iterations','referenceMacroUpdate','originalKineticResidual','coupledResidual', ...
    'collisionInvariantDefect','wallMassFlux','massError'});
writetable(summary,fullfile(folder,'summary.csv'));
writetable(summary(strcmp(summary.grid,'referenceCubic'),:),fullfile(folder,'cubic_results.csv'));
writetable(strictCIS,fullfile(folder,'strict_cis_1e9.csv'));
writetable(diagnostics,fullfile(folder,'diagnostics.csv'));
provenance=struct('date','2026-11-09','cells',64,'velocityOrder',[32 32], ...
    'mainTolerance',1e-5,'strictCISTolerance',1e-9,'stopCriterion','reference', ...
    'acceleration',.01,'gap',1,'normalization','Qp=Q/(acceleration*gap^2)=100*Q', ...
    'newGrid','NonOscillatory VELM_FLAG=3, beta=3, vx/vy both on [-4,4]', ...
    'oldGrid','vx half-Legendre on [-8,8], vy Hermite', ...
    'paper','Sharipov 1999, Table 1, printed p.398, a/b=0, BGK work [36]', ...
    'stabilization','[1 delta delta delta]', ...
    'coupling','Only the mixed collision source is updated; no H lift');
save(fullfile(folder,'comparison.mat'),'summary','strictCIS','diagnostics','provenance','runs');
fid=fopen(fullfile(folder,'provenance.json'),'w');assert(fid>=0);
fileCleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(provenance,'PrettyPrint',true));
disp(summary(strcmp(summary.grid,'referenceCubic'),:));
disp(strictCIS);
fprintf('PASS all 30 solves converged; collision invariants, wall mass flux and zero rho/T/ux perturbations verified.\n');
end

function check_run(r)
assert(r.converged && r.convergence.referenceMacroUpdate<r.config.tol);
assert(isequal(r.config.velocityOrder,[32 32]) && numel(r.velocity.weight)==1024);
assert(r.config.cells==64 && r.config.acceleration==.01);
assert(r.diagnostics.collisionInvariantDefect<1e-12 && r.diagnostics.wallMassFlux<1e-12);
assert(max(abs(r.cellMoments(:,[1 2 4])),[],'all')<1e-12);
end
