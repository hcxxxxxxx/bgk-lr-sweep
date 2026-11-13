function summary=run_gsis_cavity_comparison(savedRuns)
% Historical delta=100 comparison: preserve its delta-scaled stabilization.
% This does not use the newer Gray-aligned default unit stabilization.
% Optional two saved cavity.mat files permit analysis without repeating solves.
if nargin<1;savedRuns={};end
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
baseFolder=fullfile(root,'results','cavity_delta100_cis_tol1e5','delta_100_N32x32_stretch2');
data=load(fullfile(baseFolder,'cavity.mat'),'r');cis=data.r;
folder=fullfile(root,'results','cavity_gsis_comparison_delta100');if ~exist(folder,'dir');mkdir(folder);end
diary(fullfile(folder,'comparison.log'));cleanup=onCleanup(@()diary('off'));
rows=cell(0,13);oneStep=cell(0,6);runs=cell(1,2);
for variant=1:2
    enabled=variant==2;name='reference';if enabled;name='consistent';end
    cfg=cis.config;cfg.maxIter=500;cfg.reportEvery=20;
    cfg.macroParam=struct('kineticConsistency',enabled,'stabilization',[1 100 100 100]);
    if enabled;cfg.macroParam.relaxation=.1;end
    cfg.outputDir=fullfile(folder,name);
    [m,h]=gsis_macro_from_kinetic_2d2v(cis,cfg.macroParam);
    oneStep(end+1,:)={name,100*norm(m.nodal(:,2:3)-cis.nodalMoments(:,2:3),'fro')/norm(cis.nodalMoments(:,2:3),'fro'), ...
        h.localResidual,h.traceResidual,h.wallFluxResidual,abs(h.densityGaugeMultiplier)}; %#ok<AGROW>
    if isempty(savedRuns)
        g=monoGSISCavity2D2V(cfg);
    else
        saved=load(savedRuns{variant},'r');g=saved.r;g.reusedFromFile=savedRuns{variant};
        assert(g.config.delta==cfg.delta && isequal(g.config.cells,cfg.cells) ...
            && g.config.meshStretch==cfg.meshStretch && g.config.tol==cfg.tol ...
            && isequal(g.config.velocityOrder,cfg.velocityOrder));
        if ~isfield(g.synthetic,'relaxation')
            assert(~enabled && max(abs(g.synthetic.sourceMoments-g.synthetic.macro.nodal),[],'all')<1e-12);
            g.synthetic.relaxation=1;
        end
        assert(abs(g.synthetic.relaxation-(1-.9*enabled))<1e-14);
    end
    runs{variant}=g;
    assert(g.converged && all(isfinite(g.H),'all'));
    assert(abs(g.diagnostics.firstResidual-cis.diagnostics.firstResidual)<1e-12);
    assert(isequal(g.synthetic.hdg.stabilization,[1 100 100 100]));
    assert(g.synthetic.hdg.internalResidualCorrection==enabled);
    output=fullfile(cfg.outputDir,'delta_100_N32x32_stretch2_GSIS');
    if ~exist(output,'dir');mkdir(output);end
    if ~isempty(savedRuns)
        r=g;save(fullfile(output,'cavity.mat'),'r','-v7.3');
        writetable(array2table([g.mesh.center,g.cellMoments,g.shearStress], ...
            'VariableNames',{'x','y','rhoPert','ux','uy','TPert','Pxy'}),fullfile(output,'profile.csv'));
        writetable(array2table(g.history,'VariableNames',g.historyColumns),fullfile(output,'history.csv'));
        plotCavity2D2V(g,output);
    end
    corners=plotCavityCorners2D2V(g,output);
    candidate=corners.ellipticVelocityRoots;
    bottom=candidate(:,2)<.2 & (candidate(:,1)<.2 | candidate(:,1)>.8);
    velocityDifference=100*norm(g.cellMoments(:,2:3)-cis.cellMoments(:,2:3),'fro')/norm(cis.cellMoments(:,2:3),'fro');
    macroGap=100*norm(g.synthetic.macro.cellAverage(:,2:3)-g.cellMoments(:,2:3),'fro')/norm(g.cellMoments(:,2:3),'fro');
    rows(end+1,:)={name,cis.iterations,g.iterations,g.diagnostics.relativeKinetic,g.diagnostics.relativeCoupled, ...
        g.diagnostics.kinetic,g.diagnostics.coupled,velocityDifference,macroGap, ...
        sum(bottom),g.diagnostics.originalConverged,corners.maximumEdgeFluxProjectionError, ...
        g.synthetic.hdg.wallFluxResidual}; %#ok<AGROW>
    if enabled
        assert(g.diagnostics.originalConverged && sum(bottom)==2);
        assert(velocityDifference<.1);
    end
end
summary=cell2table(rows,'VariableNames',{'variant','cisIterations','gsisIterations', ...
    'relativeOriginalResidual','relativeCoupledResidual','absoluteOriginalResidual', ...
    'absoluteCoupledResidual','velocityDifferencePercent','macroKineticVelocityGapPercent', ...
    'bottomCornerCandidates','originalEquationPassed','streamfunctionProjectionError','macroWallResidual'});
oneStep=cell2table(oneStep,'VariableNames',{'variant','oneMacroVelocityDifferencePercent', ...
    'macroLocalResidual','macroTraceResidual','macroWallResidual','densityGaugeMultiplier'});
writetable(summary,fullfile(folder,'summary.csv'));writetable(oneStep,fullfile(folder,'one_macro.csv'));
save(fullfile(folder,'comparison.mat'),'runs','summary','oneStep','-v7.3');disp(summary);disp(oneStep);
end
