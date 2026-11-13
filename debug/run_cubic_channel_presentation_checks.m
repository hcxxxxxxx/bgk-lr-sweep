function cases = run_cubic_channel_presentation_checks
% Retest every current report metric with the reference cubic 32x32 grid.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
folder=fullfile(root,'results','channel_cubic_presentation_20261109');
if ~exist(folder,'dir');mkdir(folder);end
base=struct('deltas',[.1 1 10],'velocityOrder',[32 32], ...
    'normalQuadrature','referenceCubic','normalCutoff',4, ...
    'tangentialQuadrature','referenceCubic','tangentialCutoff',4, ...
    'reference',false,'stopCriterion','reference','tol',1e-5,'reportEvery',1000);
cases=struct('flow',{},'delta',{},'cells',{},'cis',{},'gsis',{});
refinement=table;
for N=[32 64 128]
    target=fullfile(folder,sprintf('N%d',N));
    opt=base;opt.cells=N;opt.outputDir=target;
    if isfile(fullfile(target,'comparison.mat'))
        d=load(fullfile(target,'comparison.mat'),'summary','runs');
        summary=d.summary;runs=d.runs;
    else
        [summary,runs]=runGSISChannelComparison(opt);
    end
    for i=1:numel(runs)
        for flow={'couette','poiseuille'}
            c=runs(i).cis.(flow{1});g=runs(i).gsis.(flow{1});
            check_case(c,g,N,[32 32]);
            cases(end+1)=struct('flow',flow{1},'delta',c.config.delta,'cells',N,'cis',c,'gsis',g); %#ok<AGROW>
        end
    end
    refinement=[refinement;summary]; %#ok<AGROW>
end
refinement=sortrows(refinement,{'flow','delta','cells'});
provenance=struct('date','2026-11-09','grid','Both velocities: reference cubic, beta=3, cutoff=4, 32 points', ...
    'stop','Reference-style relative macro update 1e-5', ...
    'stabilization','[1 delta delta delta]', ...
    'feedback','Qsource=Qkin+omega*(Qmacro-Qkin), no H lift', ...
    'mixing','omega=min(1,2*delta*dx/sqrt(pi))', ...
    'comparison','Old half-Legendre/Hermite results are retained separately', ...
    'independentReference','Not computed or plotted in this retest');
save(fullfile(folder,'presentation_data.mat'),'cases','refinement','provenance');
writetable(refinement,fullfile(folder,'spatial_refinement.csv'));

% Velocity quadrature sensitivity: only Nvx, then only Nvy, at fixed Nx64.
velocityCases=cases([cases.cells]==64);
for order=[64 32;128 32;128 64]'
    target=fullfile(folder,sprintf('V%d_%d',order(1),order(2)));
    opt=base;opt.cells=64;opt.velocityOrder=order';opt.outputDir=target;
    if isfile(fullfile(target,'comparison.mat'))
        d=load(fullfile(target,'comparison.mat'),'runs');runs=d.runs;
    else
        [~,runs]=runGSISChannelComparison(opt);
    end
    for i=1:numel(runs)
        for flow={'couette','poiseuille'}
            c=runs(i).cis.(flow{1});g=runs(i).gsis.(flow{1});
            check_case(c,g,64,order');
            velocityCases(end+1)=struct('flow',flow{1},'delta',c.config.delta,'cells',64,'cis',c,'gsis',g); %#ok<AGROW>
        end
    end
end
rows=cell(0,13);
for i=1:numel(velocityCases)
    vcase=velocityCases(i);c=vcase.cis;g=vcase.gsis;order=c.config.velocityOrder;
    ids=find([velocityCases.delta]==vcase.delta & strcmp({velocityCases.flow},vcase.flow));
    fine=[];
    for j=ids
        if isequal(velocityCases(j).cis.config.velocityOrder,[128 64]);fine=velocityCases(j).cis;break;end
    end
    assert(~isempty(fine));
    flowChange=NaN;qp=[NaN NaN];
    if strcmp(vcase.flow,'poiseuille')
        flowChange=100*abs(c.flowRate-fine.flowRate)/abs(fine.flowRate);
        qp=[c.flowRate g.flowRate]/(c.config.acceleration*c.config.gap^2);
    end
    rows(end+1,:)={vcase.flow,vcase.delta,order(1),order(2),c.iterations,g.iterations, ...
        100*norm(c.cellMoments(:,3)-fine.cellMoments(:,3))/norm(fine.cellMoments(:,3)), ...
        flowChange,qp(1),qp(2),100*norm(g.cellMoments(:,3)-c.cellMoments(:,3))/norm(c.cellMoments(:,3)), ...
        c.diagnostics.kinetic,g.diagnostics.kinetic}; %#ok<AGROW>
end
velocitySummary=cell2table(rows,'VariableNames',{'flow','delta','Nvx','Nvy','cisIterations','gsisIterations', ...
    'cisVelocityVs128x64Percent','cisFlowVs128x64Percent','cisQp','gsisQp', ...
    'gsisVsCISVelocityPercent','cisOriginalResidual','gsisOriginalResidual'});
writetable(velocitySummary,fullfile(folder,'velocity_refinement.csv'));
save(fullfile(folder,'velocity_refinement.mat'),'velocityCases','velocitySummary');
disp(refinement(:,{'flow','delta','cells','cisIterations','gsisIterations','relativeCellUyDifference'}));
disp(velocitySummary);
fprintf('PASS all 18 spatial cases and 24 velocity cases, CIS and GSIS, converged with unchanged stabilization.\n');
end

function check_case(c,g,N,order)
assert(c.converged && g.converged);
for r={c,g}
    run=r{1};cfg=run.config;
    assert(cfg.cells==N && isequal(cfg.velocityOrder,order));
    assert(strcmp(cfg.normalQuadrature,'referenceCubic') && strcmp(cfg.tangentialQuadrature,'referenceCubic'));
    assert(cfg.normalCutoff==4 && cfg.tangentialCutoff==4 && cfg.tol==1e-5);
    assert(strcmp(cfg.stopCriterion,'reference') && cfg.initial==0 && run.convergence.referenceMacroUpdate<1e-5);
    assert(run.diagnostics.collisionInvariantDefect<1e-12 && run.diagnostics.wallMassFlux<1e-12);
    assert(max(abs(run.cellMoments(:,[1 2 4])),[],'all')<1e-12);
end
assert(isequal(g.synthetic.hdg.cache.stabilization,[1 c.config.delta*[1 1 1]]));
end
