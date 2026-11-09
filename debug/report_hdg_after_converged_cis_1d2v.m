function summary = report_hdg_after_converged_cis_1d2v(options)
%REPORT_HDG_AFTER_CONVERGED_CIS_1D2V Converged CIS -> ONE HDG -> raw field comparison.
% No GSIS outer iteration, no omega mixing, and no modification of kinetic H.
% report_hdg_after_converged_cis_1d2v(struct('deltas',[1 10],'cells',64));
if nargin<1; options=struct; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
deltas=[1 10];
if isfield(options,'deltas'); deltas=options.deltas; options=rmfield(options,'deltas'); end
folder=fullfile(root,'results','gsis_macro_1d2v','converged_cis_one_step');
if ~isfield(options,'stopCriterion') || strcmpi(options.stopCriterion,'reference')
    folder=fullfile(folder,'reference_stopping');
end
if isfield(options,'outputDir'); folder=options.outputDir; options=rmfield(options,'outputDir'); end
if ~exist(folder,'dir'); mkdir(folder); end
diary(fullfile(folder,'run.log')); cleanup=onCleanup(@() diary('off'));
options.iterationMethod='CIS'; options.plot=false; options.save=false;
if ~isfield(options,'reference'); options.reference=false; end
rows=cell(0,13); runs=struct([]);
for i=1:numel(deltas)
    options.delta=deltas(i); kinetic=mono1D2V(options);
    if isfield(kinetic,'couette'); names={'couette','poiseuille'};
    else; names={lower(kinetic.config.flow)}; end
    runs(i).delta=deltas(i);
    for j=1:numel(names)
        name=names{j};
        if isfield(kinetic,'couette'); k=kinetic.(name); else; k=kinetic; end
        assert(k.converged,'CIS must satisfy its selected stopping criterion first.');
        if strcmpi(k.config.stopCriterion,'equation')
            assert(k.diagnostics.kinetic<k.config.tol,'CIS equation residual is too large.');
        end
        [macro,aux,hdg,hot]=gsis_macro_from_kinetic_1d2v(k,k.config.macroParam);
        difference=max(abs(macro.nodal-k.nodalMoments),[],1);
        cellUy=norm(macro.cellAverage(:,3)-k.cellMoments(:,3))/norm(k.cellMoments(:,3));
        nodalUy=norm(macro.uy-k.nodalMoments(:,3))/norm(k.nodalMoments(:,3));
        macroFlow=diff(k.mesh.edges)'*macro.cellAverage(:,3);
        rows(end+1,:)={name,k.config.delta,k.config.cells,k.iterations, ...
            cellUy,nodalUy,difference(1),difference(2),difference(4), ...
            k.flowRate,macroFlow,k.diagnostics.kinetic,mat2str(hdg.cache.stabilization)}; %#ok<AGROW>
        runs(i).(name)=struct('kinetic',k,'macro',macro,'aux',aux,'hdg',hdg,'hot',hot);
    end
end
summary=cell2table(rows,'VariableNames',{'flow','delta','cells','cisIterations', ...
    'rawRelativeCellUyDifference','rawRelativeNodalUyDifference', ...
    'maximumRhoDifference','maximumUxDifference','maximumTDifference', ...
    'cisFlowRate','rawMacroFlowRate','cisResidual','stabilization'});
writetable(summary,fullfile(folder,'summary.csv'));
save(fullfile(folder,'one_step.mat'),'summary','runs');
disp(summary);
fprintf('Raw one-step HDG comparison (NO omega mixing) saved to %s\n',folder);
fprintf('Stopping: %s, tolerance %.3e; equation residuals remain separately reported.\n', ...
    k.config.stopCriterion,k.config.tol);
end
