function summary=check_gray_gsis_coupling_1d2v(kineticConsistency,stabilization)
% Paired, fresh CIS/GSIS runs; same stopping metric, mesh and velocity rule.
% check_gray_gsis_coupling_1d2v(false) isolates the Gray-adapted wall/HOT
% changes without the additional seven-moment residual load correction.
% A second argument [1 1 1 1] compares constant penalties on the same rows.
if nargin<1;kineticConsistency=true;end
if nargin<2;stabilization=[];end
validateattributes(kineticConsistency,{'logical','numeric'},{'scalar','binary'});
if ~isempty(stabilization)
    validateattributes(stabilization,{'numeric'},{'vector','numel',4,'real','finite','positive'});
    stabilization=stabilization(:)';
end
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
folder=fullfile(root,'results','gray_consistent_coupling32');
if ~kineticConsistency;folder=fullfile(root,'results','gray_uncorrected_coupling32');end
if ~isempty(stabilization)
    label=sprintf('_S%g_%g_%g_%g',stabilization);
    folder=[folder,label];
end
if ~exist(folder,'dir');mkdir(folder);end
diary(fullfile(folder,'run.log'));cleanup=onCleanup(@()diary('off'));
base=struct('flow','both','cells',64,'gap',1,'initial',0,'velocityOrder',[32 32], ...
    'normalQuadrature','referenceCubic','normalCutoff',4, ...
    'tangentialQuadrature','referenceCubic','tangentialCutoff',4, ...
    'reference',false,'plot',false,'save',false,'stopCriterion','reference', ...
    'maxIter',10000,'reportEvery',1000, ...
    'macroParam',struct('kineticConsistency',logical(kineticConsistency)));
if ~isempty(stabilization);base.macroParam.stabilization=stabilization;end
rows=cell(0,23);oneStep=cell(0,5);
runs=struct('delta',{},'tolerance',{},'cis',{},'gsis',{});
for delta=[.1 1 10]
    for tolerance=[1e-5 1e-9 1e-11]
        cfg=base;cfg.delta=delta;cfg.tol=tolerance;
        cis=mono1D2V(cfg);gsis=monoGSIS1D2V(cfg);
        runs(end+1)=struct('delta',delta,'tolerance',tolerance,'cis',cis,'gsis',gsis); %#ok<AGROW>
        for flow={'couette','poiseuille'}
            c=cis.(flow{1});g=gsis.(flow{1});
            assert(c.converged && g.converged);
            first=c.history(1,3);assert(abs(first-g.history(1,3))<1e-15);
            expectedStabilization=[1 delta delta delta];
            if ~isempty(stabilization);expectedStabilization=stabilization;end
            assert(isequal(g.synthetic.hdg.cache.stabilization,expectedStabilization));
            assert(all(g.synthetic.omega==1));
            assert(g.synthetic.hot.consistency.enabled==logical(kineticConsistency));
            if ~kineticConsistency
                assert(all(g.synthetic.hot.consistency.loadCorrection==0,'all'));
                actualLoad=zeros(size(g.synthetic.hdg.localRhs));
                for i=1:g.mesh.count
                    ids=2*i-1:2*i;force=zeros(2,7);force(:,3)=2*g.config.acceleration;
                    force(:,5:7)=g.synthetic.hot.derivative(ids,:);
                    actualLoad(:,i)=reshape(g.synthetic.hdg.cache.M(:,:,i)*force,[],1);
                end
                assert(max(abs(actualLoad-g.synthetic.hdg.localRhs),[],'all')<1e-14);
            end
            scale=norm(c.cellMoments(:,3));
            macro=g.synthetic.macro.cellAverage;
            err=100*norm(g.cellMoments(:,3)-c.cellMoments(:,3))/scale;
            rows(end+1,:)={flow{1},delta,tolerance,c.iterations,g.iterations, ...
                c.diagnostics.kinetic,g.diagnostics.kinetic,g.diagnostics.coupled, ...
                c.diagnostics.kinetic/first,g.diagnostics.kinetic/first,g.diagnostics.coupled/first, ...
                err,100*norm(macro(:,3)-c.cellMoments(:,3))/scale, ...
                100*norm(macro(:,3)-g.cellMoments(:,3))/scale,c.flowRate,g.flowRate, ...
                max(abs(g.cellMoments(:,1)-c.cellMoments(:,1))), ...
                max(abs(g.cellMoments(:,2)-c.cellMoments(:,2))), ...
                max(abs(g.cellMoments(:,4)-c.cellMoments(:,4))), ...
                g.synthetic.hdg.boundaryFluxResidual,abs(g.synthetic.hdg.bc.gaugeMultiplier), ...
                c.convergence.referenceMacroUpdate,g.convergence.referenceMacroUpdate}; %#ok<AGROW>
            assert(g.synthetic.hdg.boundaryFluxResidual<1e-11);
            assert(abs(g.synthetic.hdg.bc.gaugeMultiplier)<1e-11);
            if tolerance==1e-11
                if kineticConsistency;assert(err<1e-5 && g.diagnostics.kinetic/first<1e-8);end
                [m,~,h]=gsis_macro_from_kinetic_1d2v(c,base.macroParam);
                oneStep(end+1,:)={flow{1},delta, ...
                    100*norm(m.cellAverage(:,3)-c.cellMoments(:,3))/scale, ...
                    h.boundaryFluxResidual,h.traceResidual}; %#ok<AGROW>
            end
            fprintf('PAIR %s delta=%g tol=%.0e: steps %d/%d, original %.4e/%.4e, velocity diff %.6g%%\n', ...
                flow{1},delta,tolerance,c.iterations,g.iterations, ...
                c.diagnostics.kinetic,g.diagnostics.kinetic,err);
        end
    end
end
summary=cell2table(rows,'VariableNames',{'flow','delta','tolerance','cisIterations','gsisIterations', ...
    'cisAbsoluteOriginalResidual','gsisAbsoluteOriginalResidual','gsisAbsoluteCoupledResidual', ...
    'cisRelativeOriginalResidual','gsisRelativeOriginalResidual','gsisRelativeCoupledResidual', ...
    'kineticVelocityDiffPercent','macroVelocityDiffPercent','macroKineticVelocityGapPercent', ...
    'cisQ','gsisQ','densityMaxDiff','normalVelocityMaxDiff','temperatureMaxDiff', ...
    'wallFluxResidual','gaugeMultiplier','cisRelativeUpdate','gsisRelativeUpdate'});
oneStep=cell2table(oneStep,'VariableNames',{'flow','delta','velocityDiffPercent','wallFluxResidual','traceResidual'});
writetable(summary,fullfile(folder,'summary.csv'));
writetable(oneStep,fullfile(folder,'converged_cis_one_macro.csv'));
save(fullfile(folder,'comparison.mat'),'summary','oneStep','runs','base','-v7.3');
disp(summary(:,[1:5 9:14]));disp(oneStep);
if kineticConsistency
    fprintf('PASS paired coupled runs, unchanged stabilization, and strict fixed-point agreement.\n');
else
    fprintf('PASS paired coupled runs with zero added consistency load and unchanged stabilization.\n');
end
end
