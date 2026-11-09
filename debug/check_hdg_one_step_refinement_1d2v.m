function refinement = check_hdg_one_step_refinement_1d2v
%CHECK_HDG_ONE_STEP_REFINEMENT_1D2V Raw CIS -> one HDG, without source mixing.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
folder=fullfile(root,'results','gsis_macro_1d2v','converged_cis_one_step');
refinement=readtable(fullfile(folder,'summary.csv'));
for N=[32 128 256]
    caseFolder=fullfile(folder,sprintf('N%d',N));
    current=report_hdg_after_converged_cis_1d2v(struct('cells',N,'outputDir',caseFolder, ...
        'stopCriterion','equation','tol',1e-10));
    audit_hdg_bulk_interface_1d2v(caseFolder);
    refinement=[refinement;current]; %#ok<AGROW>
end
refinement=sortrows(refinement,{'flow','delta','cells'});
refinement.observedOrder=nan(height(refinement),1);
for i=2:height(refinement)
    if strcmp(refinement.flow{i},refinement.flow{i-1}) && refinement.delta(i)==refinement.delta(i-1)
        refinement.observedOrder(i)=log(refinement.rawRelativeCellUyDifference(i-1) ...
            /refinement.rawRelativeCellUyDifference(i))/log(refinement.cells(i)/refinement.cells(i-1));
    end
end
writetable(refinement,fullfile(folder,'one_step_refinement.csv'));
disp(refinement(:,{'flow','delta','cells','rawRelativeCellUyDifference','observedOrder'}));
end
