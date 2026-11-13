function summary=run_gray_aligned_cavity_2d2v
% Fresh results for the latest Gray-aligned BGK algorithm; no old data reuse.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
folder=fullfile(root,'results','cavity_gsis_gray_defaults_2d2v');
if ~exist(folder,'dir');mkdir(folder);end
diary(fullfile(folder,'run.log'));cleanup=onCleanup(@()diary('off'));
rows=cell(0,11);
for delta=[100 1000]
    cfg=struct('cells',32,'meshStretch',2,'delta',delta,'velocityOrder',[32 32], ...
        'tol',1e-5,'maxIter',300,'reportEvery',20,'plot',false,'outputDir',folder);
    r=monoGSISCavity2D2V(cfg);
    assert(r.converged && ~r.synthetic.hdg.internalResidualCorrection);
    assert(r.synthetic.relaxation==1 && isequal(r.synthetic.hdg.stabilization,ones(1,4)));
    assert(strcmp(r.convergenceCriterion,'referenceCellIntegratedMacroUpdate'));
    assert(max(abs(r.synthetic.sourceMoments-r.synthetic.macro.nodal),[],'all')<1e-12);
    assert(max([r.synthetic.hdg.localResidual,r.synthetic.hdg.traceResidual, ...
        r.synthetic.hdg.wallFluxResidual])<1e-10);
    difference=NaN;
    if delta==100
        path=fullfile(root,'results','cavity_delta100_cis_tol1e5','delta_100_N32x32_stretch2','cavity.mat');
        if exist(path,'file')
            saved=load(path,'r');cis=saved.r;
            assert(isequal(cis.mesh.nodes,r.mesh.nodes) && isequal(cis.velocity.vx,r.velocity.vx));
            difference=100*norm(r.cellMoments(:,2:3)-cis.cellMoments(:,2:3),'fro')/norm(cis.cellMoments(:,2:3),'fro');
        end
    end
    gap=100*norm(r.synthetic.macro.cellAverage(:,2:3)-r.cellMoments(:,2:3),'fro')/norm(r.cellMoments(:,2:3),'fro');
    rows(end+1,:)={delta,r.iterations,r.diagnostics.macroUpdate,r.diagnostics.relativeKinetic, ...
        r.diagnostics.relativeCoupled,r.synthetic.hdg.wallFluxResidual,gap,difference, ...
        r.diagnostics.originalConverged,r.wallTime,mat2str(r.synthetic.hdg.stabilization)}; %#ok<AGROW>
end
summary=cell2table(rows,'VariableNames',{'delta','iterations','referenceMacroUpdate', ...
    'relativeOriginalResidual','relativeCoupledResidual','macroWallResidual', ...
    'macroKineticVelocityGapPercent','velocityDifferenceVsCISPercent', ...
    'originalEquationPassed','seconds','stabilization'});
writetable(summary,fullfile(folder,'summary.csv'));disp(summary);
end
