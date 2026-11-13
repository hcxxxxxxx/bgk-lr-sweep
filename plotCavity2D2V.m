function plotCavity2D2V(r,folder)
% Cell-centre streamlines/quiver and centreline profiles; no no-slip overwrite.
if nargin<2;folder='';end
nx=r.mesh.cells(1);ny=r.mesh.cells(2);U=r.config.lidSpeed;
if ~isfield(r.mesh,'xEdges')
    r.mesh.xEdges=linspace(0,1,nx+1);r.mesh.yEdges=linspace(0,1,ny+1);
end
x=(r.mesh.xEdges(1:end-1)+r.mesh.xEdges(2:end))/2;
y=(r.mesh.yEdges(1:end-1)+r.mesh.yEdges(2:end))/2;
ux=reshape(r.cellMoments(:,2),nx,ny)';uy=reshape(r.cellMoments(:,3),nx,ny)';
[X,Y]=meshgrid(x,y);speed=hypot(ux,uy);scale=max(abs(U),eps);
f=figure('Color','w','Visible','off','Position',[100 100 1600 450]);
tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
nexttile;contourf(X,Y,speed/scale,24,'LineColor','none');hold on;colorbar;
stride=max(1,floor(nx/12));ids=1:stride:nx;jds=1:stride:ny;
quiver(X(jds,ids),Y(jds,ids),ux(jds,ids)/scale,uy(jds,ids)/scale,'k');axis equal tight;
xlim([0 1]);ylim([0 1]);xlabel('x');ylabel('y');title('|u| / |U_{lid}| and velocity');
nexttile;streamslice(X,Y,ux,uy,1.5);axis equal;xlim([0 1]);ylim([0 1]);
hold on;quiver(.2,1,.6,0,0,'r','LineWidth',2,'MaxHeadSize',.3);
xlabel('x');ylabel('y');title('Cell-centre streamlines; top lid moves right');
nexttile;uxMid=interp1(x,ux',.5)';uyMid=interp1(y,uy,.5);
plot(uxMid/scale,y,'LineWidth',1.5);grid on;xlabel('u_x / |U_{lid}|');ylabel('y');
title('Vertical centreline x=0.5');ylim([0 1]);
nexttile;plot(x,uyMid/scale,'LineWidth',1.5);grid on;xlabel('x');ylabel('u_y / |U_{lid}|');
title('Horizontal centreline y=0.5');xlim([0 1]);
method='CIS';if isfield(r.config,'iterationMethod');method=r.config.iterationMethod;end
sgtitle(sprintf('%s linear BGK Q1-DG cavity: delta=%g, %dx%d cells, %dx%d velocities, %d sweeps', ...
    method,r.config.delta,r.mesh.cells,r.config.velocityOrder,r.iterations));
if ~isempty(folder)
    exportgraphics(f,fullfile(folder,'cavity.png'),'Resolution',180);
    writetable(table(y',uxMid/scale,'VariableNames',{'y','uxOverLid'}),fullfile(folder,'vertical_centerline.csv'));
    writetable(table(x',uyMid'/scale,'VariableNames',{'x','uyOverLid'}),fullfile(folder,'horizontal_centerline.csv'));
else
    set(f,'Visible','on');
end
if ~isempty(folder);close(f);end
end
