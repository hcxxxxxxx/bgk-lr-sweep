function summary=analyze_gsis_delta1000
% Reference-structure GSIS only: no internal correction, no return relaxation.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
folder=fullfile(root,'results','cavity_gsis_delta1000_reference');rows=cell(0,10);
for n=[32 64]
    stretch=2;if n==64;stretch=3;end
    path=fullfile(folder,sprintf('delta_1000_N%dx%d_stretch%g_GSIS',n,n,stretch));
    data=load(fullfile(path,'cavity.mat'),'r');r=data.r;
    assert(strcmp(r.config.iterationMethod,'GSIS') && r.config.delta==1000);
    assert(~r.synthetic.hdg.internalResidualCorrection && r.synthetic.relaxation==1);
    assert(r.converged);
    report=plotCavityCorners2D2V(r,path);
    candidate=report.ellipticVelocityRoots;
    mask=candidate(:,2)<.2 & (candidate(:,1)<.2 | candidate(:,1)>.8);
    rows(end+1,:)={n,stretch,min([diff(r.mesh.xEdges),diff(r.mesh.yEdges)]),r.iterations, ...
        r.diagnostics.relativeKinetic,r.diagnostics.relativeCoupled,sum(mask), ...
        report.maximumEdgeFluxProjectionError,mat2str(r.synthetic.hdg.stabilization),r.diagnostics.originalConverged}; %#ok<AGROW>
    zoom(r,candidate,path);
end
summary=cell2table(rows,'VariableNames',{'cellsPerDirection','meshStretch','minimumCellWidth', ...
    'iterations','relativeOriginalResidual','relativeCoupledResidual', ...
    'bottomCornerCandidates','streamfunctionProjectionError','stabilization','originalEquationPassed'});
writetable(summary,fullfile(folder,'summary.csv'));disp(summary);
end

function zoom(r,roots,folder)
nx=r.mesh.cells(1);ny=r.mesh.cells(2);
x=(r.mesh.xEdges(1:end-1)+r.mesh.xEdges(2:end))/2;
y=(r.mesh.yEdges(1:end-1)+r.mesh.yEdges(2:end))/2;
[X,Y]=meshgrid(x,y);U=reshape(r.cellMoments(:,2),nx,ny)';V=reshape(r.cellMoments(:,3),nx,ny)';
f=figure('Visible','off','Color','w','Position',[100 100 1100 900]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for row=1:2
    extent=.08;if row==2;extent=.012;end
    for side=1:2
        nexttile;
        distance=logspace(log10(x(1)*1.01),log10(extent),10);
        distanceY=logspace(log10(y(1)*1.01),log10(extent),10);
        sx=distance;bounds=[0 extent];if side==2;sx=1-distance;bounds=[1-extent 1];end
        [sx,sy]=meshgrid(sx,distanceY);
        curves=[stream2(X,Y,U,V,sx,sy,[.035 5000]);stream2(X,Y,-U,-V,sx,sy,[.035 5000])];
        lines=streamline(curves);set(lines,'Color',[0 .35 .7]);hold on;
        mask=X>=bounds(1) & X<=bounds(2) & Y<=extent;
        quiver(X(mask),Y(mask),U(mask),V(mask),'k');
        mark=roots(:,1)>bounds(1) & roots(:,1)<bounds(2) & roots(:,2)<extent;
        plot(roots(mark,1),roots(mark,2),'ro','LineWidth',1.3,'MarkerSize',7);
        axis equal;xlim(bounds);ylim([0 extent]);xlabel('x');ylabel('y');
        title(sprintf('Bottom corner %d, width=%g',side,extent));
    end
end
sgtitle(sprintf('GSIS delta=1000, %dx%d, raw velocity; no residual correction, omega=1',nx,ny));
exportgraphics(f,fullfile(folder,'corner_multiscale.png'),'Resolution',180);close(f);
end
