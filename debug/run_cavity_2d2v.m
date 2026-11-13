function summary=run_cavity_2d2v
% Fresh 24/32 cell-square calculations with a fixed 32x32 velocity rule.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
folder=fullfile(root,'results','cavity_2d2v');if ~exist(folder,'dir');mkdir(folder);end
diary(fullfile(folder,'validation.log'));cleanup=onCleanup(@()diary('off'));
rows=cell(0,11);resolution=cell(0,4);
for delta=[1 10]
    for cells=[24 32]
        r=monoCavity2D2V(struct('delta',delta,'cells',cells));
        assert(r.converged && r.diagnostics.currentWallMassFlux<1e-12);
        assert(r.diagnostics.collisionInvariantDefect<1e-12);
        assert(r.diagnostics.globalMomentumBalance<1e-8);
        U=reshape(r.cellMoments(:,2),cells,cells);V=reshape(r.cellMoments(:,3),cells,cells);
        symmetry=max([max(abs(U-flipud(U)),[],'all'),max(abs(V+flipud(V)),[],'all')]);
        assert(symmetry<1e-10);
        x=((1:cells)-.5)/cells;y=x;
        mid=interp1(x,U,.5);cross=find(mid(1:end-1)<0 & mid(2:end)>=0,1,'last');
        assert(~isempty(cross),'No central reverse-flow to lid-flow crossing.');
        yVortex=y(cross)-mid(cross)*(y(cross+1)-y(cross))/(mid(cross+1)-mid(cross));
        rows(end+1,:)={delta,cells,r.iterations,r.diagnostics.macroUpdate, ...
            r.diagnostics.kinetic,r.diagnostics.relativeKinetic, ...
            r.diagnostics.wallMassFlux,r.diagnostics.globalMomentumBalance, ...
            .5,yVortex,min(mid)/abs(r.config.lidSpeed)}; %#ok<AGROW>
        if cells==24
            coarse=r;coarseY=y;coarseMid=mid;
        else
            % Compare centreline velocities at coarse cell centres, excluding
            % corner discontinuities; this is not a mesh-independent claim.
            fineMid=interp1(y,mid,coarseY,'pchip');
            uxError=100*norm(coarseMid-fineMid)/norm(fineMid);
            vc=interp1(coarseY,reshape(coarse.cellMoments(:,3),24,24)',.5);
            vf=interp1(y,V',.5);vf=interp1(x,vf,coarseY,'pchip');
            vyError=100*norm(vc-vf)/norm(vf);
            resolution(end+1,:)={delta,uxError,vyError,yVortex}; %#ok<AGROW>
        end
    end
end
summary=cell2table(rows,'VariableNames',{'delta','cellsPerDirection','iterations', ...
    'macroRelativeUpdate','absoluteKineticResidual','relativeKineticResidual', ...
    'lastSweepWallMassFlux','globalMomentumBalance','vortexX','vortexY', ...
    'minimumCentrelineUxOverLid'});
resolution=cell2table(resolution,'VariableNames',{'delta','verticalUxChangePercent', ...
    'horizontalUyChangePercent','fineVortexY'});
writetable(summary,fullfile(folder,'summary.csv'));
writetable(resolution,fullfile(folder,'resolution.csv'));
disp(summary);disp(resolution);
end
