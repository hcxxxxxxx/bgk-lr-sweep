function check_cavity_2d2v
% Independent Gauss-integrated monolithic kinetic system, no source lag.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
cfg=struct('cells',2,'velocityOrder',[4 4],'quadrature','hermite','delta',2, ...
    'tol',1e-10,'plot',false,'save',false);
r=monoCavity2D2V(cfg);v=r.velocity;nx=2;ny=2;N=4;ns=16;nv=numel(v.weight);
A=zeros(ns*nv);b=zeros(ns*nv,1);Mglobal=zeros(ns);den=zeros(1,4);
normal=[-1 0;1 0;0 -1;0 1];cn=v.vx*normal(:,1)'+v.vy*normal(:,2)';
for side=1:4;den(side)=-sum(v.weight.*min(cn(:,side),0));end
gp=(1+[-1 1]/sqrt(3))/2;dx=1/nx;dy=1/ny;
M=zeros(4);Gx=M;Gy=M;F=cell(1,4);X=F;
for a=gp
    for c=gp
        [p,px,py]=shape(a,c,dx,dy);weight=dx*dy/4;
        M=M+weight*(p*p');Gx=Gx+weight*(px*p');Gy=Gy+weight*(py*p');
    end
end
for side=1:4
    F{side}=zeros(4);X{side}=zeros(4);
    for t=gp
        if side==1;p=shape(0,t,dx,dy);q=shape(1,t,dx,dy);len=dy;
        elseif side==2;p=shape(1,t,dx,dy);q=shape(0,t,dx,dy);len=dy;
        elseif side==3;p=shape(t,0,dx,dy);q=shape(t,1,dx,dy);len=dx;
        else;p=shape(t,1,dx,dy);q=shape(t,0,dx,dy);len=dx;end
        F{side}=F{side}+len/2*(p*p');X{side}=X{side}+len/2*(p*q');
    end
end
P=v.moment/v.gram*v.basis';
for cellId=1:N
    nodes=4*(cellId-1)+(1:4);Mglobal(nodes,nodes)=M;
    for j=1:nv
        rows=(j-1)*ns+nodes;
        A(rows,rows)=A(rows,rows)+r.config.delta*M-v.vx(j)*Gx-v.vy(j)*Gy;
        for k=1:nv
            cols=(k-1)*ns+nodes;A(rows,cols)=A(rows,cols)-r.config.delta*P(k,j)*M;
        end
        for side=1:4
            speed=cn(j,side);neighbor=r.mesh.neighbor(cellId,side);
            if speed>=0
                A(rows,rows)=A(rows,rows)+speed*F{side};
            elseif neighbor>0
                cols=(j-1)*ns+4*(neighbor-1)+(1:4);
                A(rows,cols)=A(rows,cols)+speed*X{side};
            else
                for k=find(cn(:,side)>0)'
                    cols=(k-1)*ns+nodes;
                    A(rows,cols)=A(rows,cols)+speed*v.weight(k)*cn(k,side)/den(side)*F{side};
                end
                if side==4;b(rows)=b(rows)-speed*2*r.config.lidSpeed*v.vx(j)*F{side}*ones(4,1);end
            end
        end
    end
end
original=A;rhs=b;mass=kron(v.weight',sum(Mglobal,1));
% Each velocity block is an upwind transport problem. Check all four
% quadrant scan orders independently against a sparse/direct solve.
scanError=zeros(2,2);rng(42);
for j=1:nv
    rows=(j-1)*ns+(1:ns);
    T=original(rows,rows)+r.config.delta*P(j,j)*Mglobal;
    load=randn(ns,1);direct=T\load;scanned=zeros(ns,1);
    xs=1:nx;ys=1:ny;if v.vx(j)<0;xs=nx:-1:1;end
    if v.vy(j)<0;ys=ny:-1:1;end
    for iy=ys
        for ix=xs
            c=ix+(iy-1)*nx;nodes=4*(c-1)+(1:4);
            upstream=T(nodes,:)*scanned;
            scanned(nodes)=T(nodes,nodes)\(load(nodes)-upstream);
        end
    end
    quadrant=[1+(v.vx(j)<0),1+(v.vy(j)<0)];
    scanError(quadrant(1),quadrant(2))=max(scanError(quadrant(1),quadrant(2)),max(abs(scanned-direct)));
end
assert(max(scanError,[],'all')<1e-10);
fprintf('PASS all four (vx,vy) quadrant scan orders vs direct transport solves, max %.3e.\n',max(scanError,[],'all'));
A(1,:)=mass;b(1)=0;exact=A\b;
assert(max(abs(original*exact-rhs))<1e-12);
assert(r.converged && max(abs(exact-r.H(:)))<1e-10);
assert(r.diagnostics.currentWallMassFlux<1e-13);
assert(r.diagnostics.collisionInvariantDefect<1e-13);
fprintf('PASS cavity vs independently quadrature-assembled monolithic DG: max H difference %.3e.\n', ...
    max(abs(exact-r.H(:))));
cfg.lidSpeed=-r.config.lidSpeed;reverse=monoCavity2D2V(cfg);
assert(max(abs(reverse.H+r.H),[],'all')<1e-12);
cfg.lidSpeed=0;zero=monoCavity2D2V(cfg);
assert(zero.converged && zero.iterations==1 && max(abs(zero.H),[],'all')==0);
ux=reshape(r.cellMoments(:,2),nx,ny);uy=reshape(r.cellMoments(:,3),nx,ny);
assert(max(abs(ux-flipud(ux)),[],'all')<1e-12);
assert(max(abs(uy+flipud(uy)),[],'all')<1e-12);
fprintf('PASS diffuse zero mass flux, collision invariants, forcing reversal, rest equilibrium and cavity symmetry.\n');
cfg.lidSpeed=.01;cfg.iterationMethod='GMRES';krylov=monoCavity2D2V(cfg);
assert(krylov.converged && max(abs(krylov.H-r.H),[],'all')<1e-10);
cfg.meshStretch=1.8;krylov=monoCavity2D2V(cfg);
cfg.iterationMethod='CIS';cis=monoCavity2D2V(cfg);
assert(krylov.converged && cis.converged && max(abs(krylov.H-cis.H),[],'all')<1e-10);
% Use more than two cells so tanh stretching produces a nonuniform grid.
cfg.cells=4;cfg.iterationMethod='GMRES';krylov=monoCavity2D2V(cfg);
cfg.iterationMethod='CIS';cis=monoCavity2D2V(cfg);
assert(krylov.converged && cis.converged && max(abs(krylov.H-cis.H),[],'all')<1e-10);
assert(abs(krylov.diagnostics.meanDensityPerturbation)<1e-12);
assert(krylov.diagnostics.globalMomentumBalance<1e-11);
fprintf('PASS GMRES vs CIS on uniform and stretched grids, weighted mass gauge and global momentum balance.\n');
end

function [p,px,py]=shape(x,y,dx,dy)
p=[(1-x)*(1-y);x*(1-y);(1-x)*y;x*y];
px=[-(1-y);1-y;-y;y]/dx;py=[-(1-x);-x;1-x;x]/dy;
end
