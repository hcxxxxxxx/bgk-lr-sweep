function summary = check_reference_cubic_velocity
% Compare 32x32 reference cubic quadrature with the saved Gauss baseline.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
folder=fullfile(root,'results','reference_cubic_velocity32');
if ~exist(folder,'dir');mkdir(folder);end
diary(fullfile(folder,'run.log'));cleanup=onCleanup(@()diary('off'));
data=load(fullfile(root,'results','channel_presentation_20261109','presentation_data.mat'),'cases');
cases=data.cases([data.cases.cells]==64);
variants={'original','cubicNormal','cubicBoth'};
rows=cell(0,28);runs=struct('flow',{},'delta',{},'variant',{},'cis',{},'gsis',{});strictRows=cell(0,6);
for i=1:numel(cases)
    old=cases(i);reference=old.reference.cellVelocity;
    for j=1:numel(variants)
        variant=variants{j};cfg=old.cis.config;
        cfg.reference=false;cfg.plot=false;cfg.save=false;cfg.reportEvery=1000;
        if j>1
            cfg.normalQuadrature='referenceCubic';cfg.normalCutoff=4;
        end
        if j==3
            cfg.tangentialQuadrature='referenceCubic';cfg.tangentialCutoff=4;
        end
        c=mono1D2V(cfg);g=monoGSIS1D2V(cfg);
        assert(c.converged && g.converged);
        if j==1
            assert(c.iterations==old.cis.iterations && g.iterations==old.gsis.iterations);
            assert(max(abs(c.H-old.cis.H),[],'all')<1e-14);
            assert(max(abs(g.H-old.gsis.H),[],'all')<1e-14);
        end
        for run={c,g}
            r=run{1};v=r.velocity;
            assert(numel(v.normal)==32 && numel(v.tangential)==32 && numel(v.weight)==1024);
            assert(all(v.weight>0) && all(diff(v.normal)>0));
            assert(max(abs(v.normal+flipud(v.normal)))<1e-13);
            assert(abs(sum(v.weight.*v.vx))<1e-13 && abs(sum(v.weight.*v.vy))<1e-13);
            assert(r.diagnostics.collisionInvariantDefect<1e-12);
            assert(r.diagnostics.wallMassFlux<1e-12);
            assert(max(abs(r.cellMoments(:,[1 2 4])),[],'all')<1e-12);
            if j==3
                % Independent analytic Maxwell moment checks in two velocities.
                assert(abs(sum(v.weight)-1)<1e-8);
                assert(abs(sum(v.weight.*v.vx.^2)-.5)<1e-8);
                assert(abs(sum(v.weight.*v.vy.^2)-.5)<1e-8);
                assert(abs(sum(v.weight.*v.vx.^4)-.75)<1e-6);
                assert(abs(sum(v.weight.*v.vy.^4)-.75)<1e-6);
            end
        end
        uy=c.cellMoments(:,3);ug=g.cellMoments(:,3);
        qp=[NaN NaN];paper=NaN;
        if strcmp(old.flow,'poiseuille')
            qp=[c.flowRate g.flowRate]/(cfg.acceleration*cfg.gap^2);
            ds=[.1 1 10];qs=[2.0327 1.5392 2.7673];paper=qs(ds==old.delta);
        end
        row={old.flow,old.delta,variant,64,32,32,cfg.tol,c.iterations,g.iterations, ...
            100*norm(uy-reference)/norm(reference),100*norm(ug-reference)/norm(reference), ...
            100*norm(ug-uy)/norm(uy),100*norm(g.shearStress-c.shearStress)/norm(c.shearStress), ...
            c.flowRate,g.flowRate,qp(1),qp(2),100*(qp(1)/paper-1),100*(qp(2)/paper-1), ...
            c.diagnostics.kinetic,g.diagnostics.kinetic,g.diagnostics.coupled, ...
            min(abs(c.velocity.normal)),sum(abs(c.velocity.normal)<.1), ...
            sum(c.velocity.weight)-1,sum(c.velocity.weight.*c.velocity.vy.^2)-.5, ...
            max(c.diagnostics.collisionInvariantDefect,g.diagnostics.collisionInvariantDefect), ...
            max(c.diagnostics.wallMassFlux,g.diagnostics.wallMassFlux)};
        rows(end+1,:)=row; %#ok<AGROW>
        runs(end+1)=struct('flow',old.flow,'delta',old.delta,'variant',variant,'cis',c,'gsis',g); %#ok<AGROW>
        fprintf('CHECK %s delta=%g %s: steps %d/%d; ref errors %.7g%% / %.7g%%; G-C %.7g%%\n', ...
            old.flow,old.delta,variant,c.iterations,g.iterations,row{10},row{11},row{12});
        if j==3
            cfg.tol=1e-9;tight=mono1D2V(cfg);assert(tight.converged);
            tightQp=NaN;
            if cfg.acceleration~=0;tightQp=tight.flowRate/(cfg.acceleration*cfg.gap^2);end
            strictRows(end+1,:)={old.flow,old.delta,tight.iterations, ...
                100*norm(tight.cellMoments(:,3)-reference)/norm(reference), ...
                100*norm(tight.cellMoments(:,3)-uy)/norm(tight.cellMoments(:,3)), ...
                tightQp}; %#ok<AGROW>
        end
    end
end
summary=cell2table(rows,'VariableNames',{'flow','delta','variant','cells','Nvx','Nvy', ...
    'tolerance','cisIterations','gsisIterations','cisVsBGKPercent','gsisVsBGKPercent', ...
    'relativeUyDifferencePercent','relativeShearDifferencePercent','cisQ','gsisQ','cisQp','gsisQp', ...
    'cisVsPaperPercent','gsisVsPaperPercent','cisOriginalResidual','gsisOriginalResidual', ...
    'gsisCoupledResidual','minAbsVx','countAbsVxBelow01','maxwellMassDefect', ...
    'maxwellVy2Defect','collisionInvariantDefect','wallMassFlux'});
strict=cell2table(strictRows,'VariableNames',{'flow','delta','cisIterations','cisVsBGKPercent', ...
    'toleranceEffectPercent','normalizedFlow'});
writetable(summary,fullfile(folder,'summary.csv'));
writetable(strict,fullfile(folder,'strict_stopping.csv'));
save(fullfile(folder,'comparison.mat'),'summary','strict','runs');
disp(summary(:,{'flow','delta','variant','cisIterations','gsisIterations', ...
    'cisVsBGKPercent','relativeUyDifferencePercent','cisVsPaperPercent'}));
fprintf('PASS original-grid regression, Gaussian moments, collision invariants, zero normal flow and diffuse wall mass flux.\n');
end
