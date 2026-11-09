function extend_channel_delta01_presentation
%EXTEND_CHANNEL_DELTA01_PRESENTATION Add delta=0.1 without rerunning old cases.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
folder=fullfile(root,'results','channel_presentation_20261109');
file=fullfile(folder,'presentation_data.mat');
data=load(file,'cases','refinement','provenance');
cases=data.cases;refinement=data.refinement;provenance=data.provenance;
if any([cases.delta]==0.1)
    assert(sum([cases.delta]==0.1)==6,'Incomplete cached delta=0.1 study.');
    fprintf('Six delta=0.1 cases already present.\n');return;
end
for N=[32 64 128]
    target=fullfile(folder,'delta0p1',sprintf('N%d',N));
    if isfile(fullfile(target,'comparison.mat'))
        saved=load(fullfile(target,'comparison.mat'),'runs');runs=saved.runs;
        summary=readtable(fullfile(target,'summary.csv'));
    else
        [summary,runs]=runGSISChannelComparison(struct('deltas',0.1,'cells',N, ...
            'velocityOrder',[32 32],'stopCriterion','reference','tol',1e-5, ...
            'reportEvery',1000,'outputDir',target));
    end
    for name={'couette','poiseuille'}
        c=runs.cis.(name{1});g=runs.gsis.(name{1});
        assert(c.converged && g.converged && c.config.delta==0.1 ...
            && g.config.cells==N && g.config.tol==1e-5 ...
            && strcmp(g.config.stopCriterion,'reference'));
        assert(isequal(g.synthetic.hdg.cache.stabilization,[1 .1 .1 .1]));
        n=numel(cases)+1;
        cases(n).flow=name{1};cases(n).delta=.1;cases(n).cells=N;
        cases(n).cis=c;cases(n).gsis=g;
        if N==64
            ref1=channel1D2VReference(c.config,1024,160);
            ref2=channel1D2VReference(c.config,2048,320);
            cell1=mean(reshape(ref1.velocity,[],N),1)';
            cell2=mean(reshape(ref2.velocity,[],N),1)';
            cases(n).reference=struct('cellVelocity',cell2, ...
                'coarseCellVelocity',cell1,'cells',2048,'velocityOrder',320, ...
                'selfRelativeL2',norm(cell2-cell1)/norm(cell2),'method',ref2.method);
        else
            cases(n).reference=[];
        end
    end
    refinement=[refinement;summary]; %#ok<AGROW>
end
refinement=sortrows(refinement,{'flow','delta','cells'});
provenance.delta01='Added 2026-11-09: same Nx32/64/128, Nv32x32, reference stop1e-5; unchanged solver and stabilization formula.';
save(file,'cases','refinement','provenance');
writetable(refinement,fullfile(folder,'spatial_refinement.csv'));
disp(refinement(refinement.delta==.1,{'flow','delta','cells', ...
    'cisIterations','gsisIterations','relativeCellUyDifference'}));
end
