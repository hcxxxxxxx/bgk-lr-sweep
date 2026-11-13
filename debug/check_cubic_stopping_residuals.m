function summary = check_cubic_stopping_residuals
% Separate tighter stopping from the original/coupled fixed-point residuals.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
folder=fullfile(root,'results','cubic_stopping_residuals32');
if ~exist(folder,'dir');mkdir(folder);end
diary(fullfile(folder,'run.log'));cleanup=onCleanup(@()diary('off'));
base=struct('flow','both','cells',64,'gap',1,'initial',0,'velocityOrder',[32 32], ...
    'normalQuadrature','referenceCubic','normalCutoff',4, ...
    'tangentialQuadrature','referenceCubic','tangentialCutoff',4, ...
    'reference',false,'plot',false,'save',false,'stopCriterion','reference','reportEvery',1000);
rows=cell(0,18);runs=struct('delta',{},'tolerance',{},'cis',{},'gsis',{});
for delta=[.1 1 10]
    denominator=[];
    for tolerance=[1e-5 1e-7 1e-9]
        cfg=base;cfg.delta=delta;cfg.tol=tolerance;
        c=mono1D2V(cfg);g=monoGSIS1D2V(cfg);
        runs(end+1)=struct('delta',delta,'tolerance',tolerance,'cis',c,'gsis',g); %#ok<AGROW>
        n=0;
        for flow={'couette','poiseuille'}
            n=n+1;cr=c.(flow{1});gr=g.(flow{1});
            assert(cr.converged && gr.converged);
            assert(max(cr.convergence.referenceMacroUpdate,gr.convergence.referenceMacroUpdate)<tolerance);
            first=cr.history(1,3);
            if tolerance==1e-5;denominator(n)=first;else;assert(abs(first-denominator(n))<1e-15);end %#ok<AGROW>
            assert(abs(gr.history(1,3)-first)<1e-15);
            assert(isequal(gr.synthetic.hdg.cache.stabilization,[1 delta delta delta]));
            v=gr.velocity;N=gr.mesh.count;M=kron(speye(N),gr.mesh.dx*[2 1;1 2]/6);
            sourceDifference=gr.synthetic.sourceMoments-gr.nodalMoments;
            gainDifference=(sourceDifference/v.gram)*v.basis';
            sourceResidual=delta*sqrt(sum(sum(gainDifference.*(M*gainDifference),1)'.*v.weight)/cfg.gap);
            original=gr.diagnostics.kinetic;coupled=gr.diagnostics.coupled;
            % Triangle inequality for R_original = R_coupled + source defect.
            assert(abs(original-sourceResidual)<=coupled+1e-12);
            rows(end+1,:)={flow{1},delta,tolerance,cr.iterations,gr.iterations, ...
                cr.diagnostics.kinetic/first,original/first,coupled/first, ...
                sourceResidual/first,100*norm(gr.cellMoments(:,3)-cr.cellMoments(:,3))/norm(cr.cellMoments(:,3)), ...
                cr.convergence.referenceMacroUpdate,gr.convergence.referenceMacroUpdate, ...
                cr.diagnostics.kinetic,original,coupled,first,cr.flowRate,gr.flowRate}; %#ok<AGROW>
            fprintf('STOP %s delta=%g tol=%.0e steps=%d/%d: rel CIS %.9g; rel GSIS original %.9g; coupled %.9g; source gap %.9g\n', ...
                flow{1},delta,tolerance,cr.iterations,gr.iterations, ...
                cr.diagnostics.kinetic/first,original/first,coupled/first,sourceResidual/first);
        end
    end
end
summary=cell2table(rows,'VariableNames',{'flow','delta','tolerance','cisIterations','gsisIterations', ...
    'cisRelativeOriginalResidual','gsisRelativeOriginalResidual','gsisRelativeCoupledResidual', ...
    'relativeSourceDefect','gsisVsCISVelocityPercent','cisRelativeUpdate','gsisRelativeUpdate', ...
    'cisAbsoluteOriginalResidual','gsisAbsoluteOriginalResidual','gsisAbsoluteCoupledResidual', ...
    'commonCISFirstResidual','cisQ','gsisQ'});
writetable(summary,fullfile(folder,'summary.csv'));
save(fullfile(folder,'comparison.mat'),'summary','runs');
disp(summary(:,1:10));
fprintf('PASS all stopping levels, common first-residual normalization and source-defect norm consistency.\n');
end
