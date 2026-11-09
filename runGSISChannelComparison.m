function [summary,runs] = runGSISChannelComparison(options)
%RUNGSISCHANNELCOMPARISON Compare coupled GSIS and CIS on identical grids.
% [summary,runs]=runGSISChannelComparison;
% runGSISChannelComparison(struct('deltas',[1 10],'cells',128));
% Other options are forwarded to mono1D2V. Solver plots/saves are suppressed;
% CSV/MAT/log results are written to outputDir, with separate residuals for
% the original CIS equation and the actual mixed-source coupled equation.
if nargin<1; options=struct; end
root=fileparts(mfilename('fullpath'));
deltas=[1 10];
if isfield(options,'deltas'); deltas=options.deltas; options=rmfield(options,'deltas'); end
validateattributes(deltas,{'numeric'},{'vector','real','finite','positive'});
folder=fullfile(root,'results','gsis_channel_comparison');
if ~isfield(options,'stopCriterion') || strcmpi(options.stopCriterion,'reference')
    folder=fullfile(folder,'reference_stopping');
end
if isfield(options,'outputDir'); folder=options.outputDir; options=rmfield(options,'outputDir'); end
if ~exist(folder,'dir'); mkdir(folder); end
diary(fullfile(folder,'run.log')); cleanup=onCleanup(@() diary('off'));
options.plot=false; options.save=false;
if ~isfield(options,'reference'); options.reference=false; end
rows=cell(0,20); runs=struct([]);
for i=1:numel(deltas)
    options.delta=deltas(i); options.iterationMethod='CIS';
    cis=mono1D2V(options); gsis=monoGSIS1D2V(options);
    runs(i).delta=deltas(i); runs(i).cis=cis; runs(i).gsis=gsis;
    if isfield(cis,'couette'); names={'couette','poiseuille'};
    else; names={lower(cis.config.flow)}; end
    for j=1:numel(names)
        name=names{j};
        if isfield(cis,'couette'); c=cis.(name); g=gsis.(name);
        else; c=cis; g=gsis; end
        assert(c.converged && g.converged,'Comparison contains an unconverged solve.');
        scale=norm(c.cellMoments(:,3));
        cellChange=norm(g.cellMoments(:,3)-c.cellMoments(:,3))/max(scale,eps);
        nodalChange=norm(g.nodalMoments(:,3)-c.nodalMoments(:,3))/max(norm(c.nodalMoments(:,3)),eps);
        flowChange=NaN;
        if abs(c.flowRate)>1e-12; flowChange=(g.flowRate-c.flowRate)/abs(c.flowRate); end
        s=g.synthetic.hdg.cache.stabilization;
        rows(end+1,:)={name,deltas(i),c.config.cells,c.iterations,g.iterations, ...
            cellChange,nodalChange,c.flowRate,g.flowRate,flowChange, ...
            c.diagnostics.kinetic,g.diagnostics.kinetic,g.diagnostics.coupled, ...
            min(g.synthetic.omega),max(g.synthetic.omega),mat2str(s), ...
            c.config.stopCriterion,c.config.tol,c.convergence.referenceMacroUpdate, ...
            g.convergence.referenceMacroUpdate}; %#ok<AGROW>
    end
end
summary=cell2table(rows,'VariableNames',{'flow','delta','cells','cisIterations','gsisIterations', ...
    'relativeCellUyDifference','relativeNodalUyDifference','cisFlowRate','gsisFlowRate', ...
    'signedRelativeFlowDifference','cisResidual','gsisOriginalCISResidual','gsisCoupledResidual', ...
    'minimumOmega','maximumOmega','stabilization','stopCriterion','stopTolerance', ...
    'cisReferenceMacroUpdate','gsisReferenceMacroUpdate'});
writetable(summary,fullfile(folder,'summary.csv'));
save(fullfile(folder,'comparison.mat'),'summary','runs');
disp(summary);
fprintf('CIS/GSIS comparison saved to %s\n',folder);
end
