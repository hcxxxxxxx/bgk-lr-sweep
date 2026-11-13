function summary=compare_su2019_bgk_couette(mode,reuse)
% Compare shear observables with JCP 378 (2019), Table 1, not its SIS algorithm.
% The full 3V/2V models differ; only the isothermal BGK shear reduction is
% used for the velocity/stress benchmark. No cross-algorithm speedup claim.
% Main runs retain the Gray-style BGK adaptation: S=1, no residual load,
% full macro source return, no H lift. Production solver defaults are untouched.
% Table data are the FIRST (Jiang & Luo) row at each inverse delta, rather
% than the slightly less accurate second (Su SIS) row. The columns are
% U1(1)/2, U1'(1/2)/2, -P12/4, where P12=2*<vx*vy*H>/wallSpeed.
if nargin<1;mode='baseline';end
if nargin<2;reuse=false;end
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
folder=fullfile(root,'results','su2019_bgk_couette');
if ~exist(folder,'dir');mkdir(folder);end
diary(fullfile(folder,[mode,'.log']));cleanup=onCleanup(@()diary('off'));
paper=[.01 .013430729 .025200983 .2796682147; ...
       .1  .072922113 .132195579 .2611624603; ...
       1   .251861340 .444228470 .1694625753; ...
       10  .441224641 .835285766 .04155607783; ...
       100 .493069780 .980081002 .004900405010];
writetable(array2table(paper,'VariableNames',{'delta','wallVelocityOver2','centerGradientOver2','minusP12Over4'}), ...
    fullfile(folder,'paper_table1.csv'));
base=struct('flow','Couette','gap',1,'wallVelocity',[-1 1],'acceleration',0, ...
    'cells',128,'velocityOrder',[32 32],'normalQuadrature','referenceCubic','normalCutoff',4, ...
    'tangentialQuadrature','referenceCubic','tangentialCutoff',4, ...
    'reference',false,'initial',0,'plot',false,'save',false,'maxIter',30000, ...
    'tol',1e-5,'stopCriterion','reference','checkEvery',1,'reportEvery',1000, ...
    'macroParam',struct('kineticConsistency',false,'stabilization',ones(1,4)));
cases=struct('delta',{},'cells',{},'velocityOrder',{},'tol',{},'method',{});
switch mode
    case 'baseline'
        for d=paper(:,1)'
            for method={'CIS','GSIS'}
                cases(end+1)=struct('delta',d,'cells',128,'velocityOrder',[32 32], ...
                    'tol',1e-5,'method',method{1}); %#ok<AGROW>
            end
        end
    case 'refinement'
        for d=paper(:,1)'
            for method={'CIS','GSIS'}
                cases(end+1)=struct('delta',d,'cells',128,'velocityOrder',[32 32], ...
                    'tol',1e-8,'method',method{1}); %#ok<AGROW>
            end
        end
        for d=paper(:,1)'
            cases(end+1)=struct('delta',d,'cells',512,'velocityOrder',[64 32], ...
                'tol',1e-8,'method','GSIS'); %#ok<AGROW>
        end
    case 'isolation'
        for method={'CIS','GSIS'}
            cases(end+1)=struct('delta',.01,'cells',128,'velocityOrder',[64 32], ...
                'tol',1e-8,'method',method{1}); %#ok<AGROW>
        end
        for d=[10 100]
            cases(end+1)=struct('delta',d,'cells',512,'velocityOrder',[32 32], ...
                'tol',1e-8,'method','GSIS'); %#ok<AGROW>
        end
    otherwise
        error('Choose baseline, refinement, or isolation.');
end
rows=cell(0,30);
names={'delta','method','cells','normalVelocityPoints','tangentialVelocityPoints','tolerance', ...
    'iterations','macroStopPassed','macroUpdate','wallVelocityOver2','paperWallVelocityOver2', ...
    'wallSignedErrorPercent','wallErrorPercent','centerGradientOver2','paperCenterGradientOver2', ...
    'gradientSignedErrorPercent','gradientErrorPercent','minusP12Over4','paperMinusP12Over4', ...
    'shearSignedErrorPercent','shearErrorPercent','absoluteOriginalResidual','relativeOriginalResidual', ...
    'absoluteCoupledResidual','relativeCoupledResidual','macroKineticVelocityGapPercent', ...
    'shearRelativeSpan','sourceWallVelocityOver2','interiorEndpointVelocityOver2','seconds'};
for test=cases
    cfg=base;cfg.delta=test.delta;cfg.cells=test.cells;cfg.velocityOrder=test.velocityOrder;cfg.tol=test.tol;
    tag=sprintf('delta_%g_%s_N%d_V%d_%d_tol%.0e',cfg.delta,test.method,cfg.cells,cfg.velocityOrder,cfg.tol);
    if reuse
        saved=load(fullfile(folder,[tag,'.mat']),'r');r=saved.r;
    elseif strcmp(test.method,'CIS');r=mono1D2V(cfg);else;r=monoGSIS1D2V(cfg);end
    reference=paper(paper(:,1)==cfg.delta,2:4);assert(numel(reference)==3);
    U=cfg.wallVelocity(2);u=r.nodalMoments(:,3)/U;
    slopes=(u(2:2:end)-u(1:2:end))/(cfg.gap/cfg.cells);
    center=mean(slopes(cfg.cells/2+[0 1]));
    % The physical wall trace combines prescribed incoming and interior
    % outgoing H. The incoming DG endpoint polynomial is only weakly tied
    % to the boundary condition, so its moment is a different diagnostic.
    physicalWall=r.kineticTrace(end,:)*r.velocity.moment(:,3)/U;
    value=[physicalWall/2,center/2,-mean(r.shearStress)/(2*U)];
    err=100*(value-reference)./reference;
    sourceWall=NaN;gap=0;
    if strcmp(test.method,'GSIS')
        assert(~r.synthetic.hot.consistency.enabled && all(r.synthetic.omega==1));
        assert(isequal(r.synthetic.hdg.cache.stabilization,ones(1,4)));
        sourceWall=r.synthetic.sourceMoments(end,3)/(2*U);
        gap=100*norm(r.synthetic.sourceMoments(:,3)-r.nodalMoments(:,3))/norm(r.nodalMoments(:,3));
    end
    first=r.history(1,3);span=(max(r.shearStress)-min(r.shearStress))/max(abs(mean(r.shearStress)),realmin);
    rows(end+1,:)={cfg.delta,test.method,cfg.cells,cfg.velocityOrder(1),cfg.velocityOrder(2),cfg.tol, ...
        r.iterations,r.converged,r.convergence.referenceMacroUpdate,value(1),reference(1),err(1),abs(err(1)), ...
        value(2),reference(2),err(2),abs(err(2)),value(3),reference(3),err(3),abs(err(3)), ...
        r.diagnostics.kinetic,r.diagnostics.kinetic/first,r.diagnostics.coupled,r.diagnostics.coupled/first, ...
        gap,span,sourceWall,u(end)/2,r.wallTime}; %#ok<AGROW>
    save(fullfile(folder,[tag,'.mat']),'r','reference','value','err','-v7.3');
    writetable(array2table([r.mesh.nodes+.5,u],'VariableNames',{'paperNormalCoordinate','U1'}),fullfile(folder,[tag,'_profile.csv']));
    summary=cell2table(rows,'VariableNames',names);
    writetable(summary,fullfile(folder,[mode,'_summary.csv']));
    fprintf('PAPER CHECK delta=%g %s Nx=%d V=%dx%d tol=%.0e: steps=%d, wall/gradient/shear error=[%.6g %.6g %.6g] percent, original rel=%.4e\n', ...
        cfg.delta,test.method,cfg.cells,cfg.velocityOrder,cfg.tol,r.iterations,err,r.diagnostics.kinetic/first);
end
disp(summary(:,{'delta','method','cells','normalVelocityPoints','tolerance','iterations','wallErrorPercent','gradientErrorPercent','shearErrorPercent','relativeOriginalResidual'}));
end
