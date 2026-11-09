function cases = run_channel_presentation_checks
%RUN_CHANNEL_PRESENTATION_CHECKS Common stopping/grid data for plotting.
% Does not change stabilization, coupling, or existing saved solver results.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
folder=fullfile(root,'results','channel_presentation_20261109');
if ~exist(folder,'dir'); mkdir(folder); end
cases=struct([]); refinement=table; n=0;
for N=[32 64 128]
    target=fullfile(folder,sprintf('N%d',N));
    if ~exist(target,'dir'); mkdir(target); end
    dataFile=fullfile(target,'comparison.mat');
    if N==64 && ~isfile(dataFile)
        baseline=fullfile(root,'results','gsis_channel_comparison','reference_stopping');
        copyfile(fullfile(baseline,'comparison.mat'),dataFile);
        copyfile(fullfile(baseline,'summary.csv'),fullfile(target,'summary.csv'));
    end
    if isfile(dataFile)
        data=load(dataFile,'runs'); runs=data.runs;
        summary=readtable(fullfile(target,'summary.csv'));
    else
        [summary,runs]=runGSISChannelComparison(struct('cells',N, ...
            'velocityOrder',[32 32],'stopCriterion','reference','tol',1e-5, ...
            'reportEvery',1000,'outputDir',target));
    end
    for i=1:numel(runs)
        for name={'couette','poiseuille'}
            c=runs(i).cis.(name{1}); g=runs(i).gsis.(name{1});
            for r={c,g}
                cfg=r{1}.config;
                assert(r{1}.converged && cfg.cells==N && cfg.tol==1e-5 ...
                    && strcmp(cfg.stopCriterion,'reference') ...
                    && isequal(cfg.velocityOrder,[32 32]));
                assert(strcmp(cfg.normalQuadrature,'halfLegendre') && cfg.initial==0);
            end
            assert(isequal(g.synthetic.hdg.cache.stabilization,[1 c.config.delta*[1 1 1]]));
            n=n+1;
            cases(n).flow=name{1}; cases(n).delta=c.config.delta; cases(n).cells=N;
            cases(n).cis=c; cases(n).gsis=g;
            if N==64
                % Independent BGK reference, with a separate self-refinement
                % check. These values are not the published BGK table.
                cfg=c.config;
                ref1=channel1D2VReference(cfg,1024,160);
                ref2=channel1D2VReference(cfg,2048,320);
                cell1=mean(reshape(ref1.velocity,[],N),1)';
                cell2=mean(reshape(ref2.velocity,[],N),1)';
                cases(n).reference=struct('cellVelocity',cell2, ...
                    'coarseCellVelocity',cell1,'cells',2048,'velocityOrder',320, ...
                    'selfRelativeL2',norm(cell2-cell1)/norm(cell2), ...
                    'method',ref2.method);
            else
                cases(n).reference=[];
            end
        end
    end
    refinement=[refinement;summary]; %#ok<AGROW>
end
refinement=sortrows(refinement,{'flow','delta','cells'});
writetable(refinement,fullfile(folder,'spatial_refinement.csv'));
provenance=struct('date','2026-11-09','source', ...
    'N64 copied from Nov 5 reference_stopping; N32/N128 computed with identical stopping', ...
    'grid','P1; Nv=32x32; half-Legendre normal [-8,8]; tangential Hermite', ...
    'stop','max relative update of cell integrals of physical rho,T,and velocity magnitude < 1e-5', ...
    'sourceMonitor','CIS kinetic moments; GSIS mixed source moments', ...
    'stabilization','[1 delta delta delta]', ...
    'mixing','omega=min(1,2*delta*dx/sqrt(pi)); changes with spatial resolution', ...
    'independentReference','Characteristic integral BGK; 2048 cells, 320 positive speeds; not literature values');
save(fullfile(folder,'presentation_data.mat'),'cases','refinement','provenance');
disp(refinement(:,{'flow','delta','cells','cisIterations','gsisIterations','relativeCellUyDifference'}));
fprintf('Presentation data saved to %s\n',folder);
end
