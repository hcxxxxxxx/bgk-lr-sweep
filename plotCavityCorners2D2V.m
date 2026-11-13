function report=plotCavityCorners2D2V(r,folder)
% Raw cell-centre velocity streamlines and signed streamfunction diagnostics.
% No artificial corner eddies or wall no-slip velocities are inserted.
if nargin<2;folder=r.config.outputDir;end
if ~exist(folder,'dir');mkdir(folder);end
nx=r.mesh.cells(1);ny=r.mesh.cells(2);v=r.velocity;H=r.H;
if ~isfield(r.mesh,'xEdges')
    r.mesh.xEdges=linspace(0,1,nx+1);r.mesh.yEdges=linspace(0,1,ny+1);
end
x=(r.mesh.xEdges(1:end-1)+r.mesh.xEdges(2:end))/2;
y=(r.mesh.yEdges(1:end-1)+r.mesh.yEdges(2:end))/2;
U=reshape(r.cellMoments(:,2),nx,ny)';V=reshape(r.cellMoments(:,3),nx,ny)';
[X,Y]=meshgrid(x,y);lid=abs(r.config.lidSpeed);
% Edge-integrated kinetic mass flux determines differences of psi at
% vertices. Least squares removes the small remaining mass residual.
% The values on all four impermeable walls are fixed to psi=0.
nvertex=(nx+1)*(ny+1);rows=[];cols=[];vals=[];target=[];edge=0;
for iy=1:ny
    for ix=2:nx
        left=ix-1+(iy-1)*nx;right=ix+(iy-1)*nx;
        trace=H(4*(left-1)+[2 4],:);incoming=v.vx<0;
        trace(:,incoming)=H(4*(right-1)+[1 3],incoming);
        flux=sum(trace*(v.weight.*v.vx))*diff(r.mesh.yEdges(iy:iy+1))/2;
        edge=edge+1;a=ix+(iy-1)*(nx+1);b=ix+iy*(nx+1);
        rows=[rows;edge;edge];cols=[cols;a;b];vals=[vals;-1;1];target(edge,1)=flux; %#ok<AGROW>
    end
end
for iy=2:ny
    for ix=1:nx
        low=ix+(iy-2)*nx;high=ix+(iy-1)*nx;
        trace=H(4*(low-1)+[3 4],:);incoming=v.vy<0;
        trace(:,incoming)=H(4*(high-1)+[1 2],incoming);
        flux=-sum(trace*(v.weight.*v.vy))*diff(r.mesh.xEdges(ix:ix+1))/2;
        edge=edge+1;a=ix+(iy-1)*(nx+1);b=a+1;
        rows=[rows;edge;edge];cols=[cols;a;b];vals=[vals;-1;1];target(edge,1)=flux; %#ok<AGROW>
    end
end
incidence=sparse(rows,cols,vals,edge,nvertex);
[ix,iy]=ndgrid(1:nx+1,1:ny+1);free=ix>1 & ix<nx+1 & iy>1 & iy<ny+1;
psi=zeros(nvertex,1);psi(free)=incidence(:,free)\target;
fluxMismatch=incidence*psi-target;
P=reshape(psi,nx+1,ny+1)';
% Also integrate actual cell velocities vertically. It is deliberately a
% separate diagnostic, since DG velocity is not pointwise divergence free.
q=reshape(r.nodalMoments(:,2),4,nx,ny);
low=reshape((q(1,:,:)+q(2,:,:))/2,nx,ny)';
high=reshape((q(3,:,:)+q(4,:,:))/2,nx,ny)';
height=diff(r.mesh.yEdges(:));increment=height.*(low+high)/2;
pathPsi=[zeros(1,nx);cumsum(increment,1)];
pathCentre=pathPsi(1:end-1,:)+height.*(3*low+high)/8;
% Find roots of BOTH velocity components in each bilinear interpolant.
% Report a candidate only if its linearized vector field is elliptic.
candidates=[];
for j=1:ny-1
    for i=1:nx-1
        uv=[U(j,i),U(j,i+1),U(j+1,i),U(j+1,i+1)];
        vv=[V(j,i),V(j,i+1),V(j+1,i),V(j+1,i+1)];
        if min(uv)>0 || max(uv)<0 || min(vv)>0 || max(vv)<0;continue;end
        a=[uv(1),uv(2)-uv(1),uv(3)-uv(1),uv(4)-uv(3)-uv(2)+uv(1)];
        b=[vv(1),vv(2)-vv(1),vv(3)-vv(1),vv(4)-vv(3)-vv(2)+vv(1)];
        z=[.5;.5];
        for step=1:15
            f=[a(1)+a(2)*z(1)+a(3)*z(2)+a(4)*prod(z);b(1)+b(2)*z(1)+b(3)*z(2)+b(4)*prod(z)];
            J=[a(2)+a(4)*z(2),a(3)+a(4)*z(1);b(2)+b(4)*z(2),b(3)+b(4)*z(1)];
            if rcond(J)<1e-12;break;end
            dz=J\f;z=z-dz;if norm(dz)<1e-12;break;end
        end
        if any(z<0) || any(z>1) || norm(f)>lid*1e-8;continue;end
        J=J*diag(1./[x(i+1)-x(i),y(j+1)-y(j)]);
        if sum(diag(J))^2-4*det(J)<0
            px=x(i)+z(1)*(x(i+1)-x(i));py=y(j)+z(2)*(y(j+1)-y(j));
            candidates(end+1,:)=[px,py,J(2,1)-J(1,2),interp2(r.mesh.xEdges,r.mesh.yEdges,P,px,py)]; %#ok<AGROW>
        end
    end
end
if isempty(candidates);candidates=zeros(0,4);end
if ~isempty(candidates)
    % Deduplicate positions only; rounding psi would erase higher-order eddies.
    [~,keep]=unique(round(candidates(:,1:2),10),'rows');candidates=candidates(keep,:);
end
f=figure('Visible','off','Color','w','Position',[100 100 1400 800]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile;streamslice(X,Y,U,V,2);axis equal;xlim([0 1]);ylim([0 1]);
title('Full raw velocity streamlines');xlabel('x');ylabel('y');
for side=1:2
    nexttile;if side==1;bounds=[0 .18];else;bounds=[.82 1];end
    % Seeds are placed in the bottom corners; integrate in both directions.
    [sx,sy]=meshgrid(linspace(max(bounds(1),x(1)),min(bounds(2),x(end)),9), ...
        linspace(y(1),.18,9));
    curves=[stream2(X,Y,U,V,sx,sy,[.05 5000]);stream2(X,Y,-U,-V,sx,sy,[.05 5000])];
    h=streamline(curves);set(h,'Color',[0 .35 .7]);hold on;
    mask=X>=bounds(1) & X<=bounds(2) & Y<=.18;
    quiver(X(mask),Y(mask),U(mask),V(mask),'k');
    marks=candidates(:,1)>=bounds(1) & candidates(:,1)<=bounds(2) & candidates(:,2)<=.18;
    plot(candidates(marks,1),candidates(marks,2),'ro','MarkerSize',7,'LineWidth',1.2);
    axis equal;xlim(bounds);ylim([0 .18]);xlabel('x');ylabel('y');
    title(sprintf('Bottom %s: raw velocity',string_side(side)));
end
nexttile;contour(r.mesh.xEdges,r.mesh.yEdges,P/lid,30);axis equal tight;colorbar;
title('Psi / |U_lid| from kinetic face fluxes');xlabel('x');ylabel('y');
for side=1:2
    nexttile;if side==1;bounds=[0 .18];else;bounds=[.82 1];end
    mask=r.mesh.xEdges>=bounds(1) & r.mesh.xEdges<=bounds(2);ymask=r.mesh.yEdges<=.18;
    range=P(ymask,mask)/lid;level=max(abs(range),[],'all');
    levels=unique([-logspace(-10,log10(max(level,1e-10)),24),0,logspace(-10,log10(max(level,1e-10)),24)]);
    contour(r.mesh.xEdges,r.mesh.yEdges,P/lid,sort(levels));hold on;
    if any(range>0,'all') && any(range<0,'all')
        contour(r.mesh.xEdges,r.mesh.yEdges,P/lid,[0 0],'k','LineWidth',1.5);
    end
    marks=candidates(:,1)>=bounds(1) & candidates(:,1)<=bounds(2) & candidates(:,2)<=.18;
    plot(candidates(marks,1),candidates(marks,2),'ro','MarkerSize',7,'LineWidth',1.2);
    axis equal;xlim(bounds);ylim([0 .18]);xlabel('x');ylabel('y');
    title(sprintf('Bottom %s: signed Psi',string_side(side)));
end
method='CIS';if isfield(r.config,'iterationMethod');method=r.config.iterationMethod;end
sgtitle(sprintf('%s cavity delta=%g, %dx%d, tol=%g; raw kinetic velocity',method,r.config.delta,nx,ny,r.config.tol));
exportgraphics(f,fullfile(folder,'bottom_corners.png'),'Resolution',180);close(f);
report=struct('ellipticVelocityRoots',candidates,'maximumEdgeFluxProjectionError',max(abs(fluxMismatch)), ...
    'relativeEdgeFluxProjectionError',norm(fluxMismatch)/max(norm(target),eps), ...
    'psi',P,'verticalPathPsi',pathCentre,'config',r.config,'diagnostics',r.diagnostics);
save(fullfile(folder,'corner_diagnostics.mat'),'report');
writetable(array2table(candidates,'VariableNames',{'x','y','vorticity','psi'}),fullfile(folder,'vortex_candidates.csv'));
disp(report.ellipticVelocityRoots);
fprintf('Psi edge-flux projection error: absolute %.4e, relative %.4e\n', ...
    report.maximumEdgeFluxProjectionError,report.relativeEdgeFluxProjectionError);
end

function label=string_side(side)
if side==1;label='left';else;label='right';end
end
