function c=hdg_precompute_q1_2d2v(mesh,v,delta,param)
% Q1 HDG for 2D2V linear BGK. Cell fields [p ux uy T Lxx Lxy Ex Ey].
% sigma_xx=-Lxx/delta, sigma_xy=-Lxy/delta, q=-E/delta.
% Each face has two trace nodes for each of [p ux uy T].
% Gray archive mapping: p corresponds to its scalar T, u to its vector q.
% Keep the cell normal velocity in the first flux and common pressure trace
% in momentum, as CELL_FLUX/TRACE_FLUX in Synthetic_Acceleration.f90.
% Stress/heat factors below are BGK 2V coefficients, not Callaway constants.
if nargin<4;param=struct;end
s=ones(1,4);if isfield(param,'stabilization');s=param.stabilization(:)';end
validateattributes(s,{'numeric'},{'numel',4,'positive','finite','real'});
nx=mesh.cells(1);ny=mesh.cells(2);N=mesh.count;nf=(nx+1)*ny+nx*(ny+1);nt=8*nf;
m=[2 1;1 2]/6;g=[-1 -1;1 1]/2;I=eye(4);
E={I(:,[1 3]),I(:,[2 4]),I(:,[1 2]),I(:,[3 4])};normals=[-1 0;1 0;0 -1;0 1];
Fx=zeros(4,8);Fy=Fx;Fx(1,2)=1;Fy(1,3)=1;
Fx(2,[1 5])=[1 -1/delta];Fy(2,6)=-1/delta;
Fx(3,6)=-1/delta;Fy(3,[1 5])=[1 1/delta];
Fx(4,7)=-1/delta;Fy(4,8)=-1/delta;
Ax=zeros(4,8);Ay=Ax;Ax(1,2)=1;Ay(1,3)=-1;
Ax(2,3)=1;Ay(2,2)=1;Ax(3,4)=1;Ay(4,4)=1;
c=struct('formulation','gray-bgk-2v-v2','count',N,'faceCount',nf,'delta',delta,'stabilization',s,'E',{E}, ...
    'normals',normals,'A',zeros(32,32,N),'B',zeros(32,32,N), ...
    'C',zeros(8,32,4,N),'D',zeros(8,8,4,N), ...
    'CJ',zeros(8,32,4,N),'response',zeros(32,32,N), ...
    'ids',zeros(32,N),'faces',zeros(N,4),'M',zeros(4,4,N), ...
    'Gx',zeros(4,4,N),'Gy',zeros(4,4,N),'faceM',zeros(2,2,4,N), ...
    'wallJ',zeros(4,4,4),'edgesX',mesh.xEdges,'edgesY',mesh.yEdges);
fluxBasis=[ones(size(v.vx)),2*v.vx,2*v.vy,v.speed2-2];
D=[1 0 0 0;0 1 0 0;0 0 1 0;-1 0 0 1];tangent=(D/v.gram*v.basis')';
for side=1:4
    cn=v.vx*normals(side,1)+v.vy*normals(side,2);out=cn>0;in=cn<0;W=tangent;
    reflected=-sum((v.weight(out).*cn(out)).*tangent(out,:),1)/sum(v.weight(in).*cn(in));
    W(in,:)=repmat(reflected,nnz(in),1);
    c.wallJ(:,:,side)=fluxBasis'*((v.weight.*cn).*W);
end
rows=zeros(1024,N);cols=rows;values=rows;meanRow=zeros(1,nt);
for j=1:N
    ix=mod(j-1,nx)+1;iy=floor((j-1)/nx)+1;
    faces=[ix+(iy-1)*(nx+1),ix+1+(iy-1)*(nx+1), ...
        (nx+1)*ny+ix+(iy-1)*nx,(nx+1)*ny+ix+iy*nx];
    ids=reshape((8*(faces-1)+(1:8)') ,[],1);c.ids(:,j)=ids;c.faces(j,:)=faces;
    dx=mesh.dx(j);dy=mesh.dy(j);M=dx*dy*kron(m,m);Gx=dy*kron(m,g);Gy=dx*kron(g,m);
    A=[-kron(Fx,Gx)-kron(Fy,Gy);kron(Ax,Gx)+kron(Ay,Gy)];
    A(17:32,17:32)=A(17:32,17:32)+kron(eye(4),M);B=zeros(32);
    for side=1:4
        normal=normals(side,:);length=dy;if side>2;length=dx;end
        Mf=length*m;F=normal(1)*Fx+normal(2)*Fy;
        tracePressure=zeros(4);tracePressure(2:3,1)=normal';F(2:3,1)=0;
        C=kron(F+[diag(s),zeros(4)],Mf*E{side}');Dface=kron(tracePressure-diag(s),Mf);
        localRows=8*(side-1)+(1:8);
        A(1:16,:)=A(1:16,:)+kron(eye(4),E{side})*C;
        B(1:16,localRows)=-kron(eye(4),E{side})*Dface;
        B(17:32,localRows)=kron(normal(1)*Ax(:,1:4)+normal(2)*Ay(:,1:4),E{side}*Mf);
        c.C(:,:,side,j)=C;c.D(:,:,side,j)=Dface;c.CJ(:,:,side,j)=C;
        if mesh.neighbor(j,side)==0
            c.CJ(:,:,side,j)=C-[kron(c.wallJ(:,:,side),Mf*E{side}'),zeros(8,16)];
        end
        c.faceM(:,:,side,j)=Mf;
    end
    response=A\B;K=zeros(32);
    for side=1:4
        localRows=8*(side-1)+(1:8);
        K(localRows,:)=c.CJ(:,:,side,j)*response;
        K(localRows,localRows)=K(localRows,localRows)+c.D(:,:,side,j);
    end
    density=zeros(1,32);density(1:4)=sum(M,1);density(13:16)=-sum(M,1);
    meanRow(ids)=meanRow(ids)+density*response;
    [rr,cc]=ndgrid(ids,ids);rows(:,j)=rr(:);cols(:,j)=cc(:);values(:,j)=K(:);
    c.A(:,:,j)=A;c.B(:,:,j)=B;c.response(:,:,j)=response;
    c.M(:,:,j)=M;c.Gx(:,:,j)=Gx;c.Gy(:,:,j)=Gy;
end
c.K=sparse(rows(:),cols(:),values(:),nt,nt);
compatibility=zeros(nt,1);compatibility(1:8:end)=1;compatibility(2:8:end)=1;
compatibility=compatibility/norm(compatibility);
c.meanRow=meanRow;c.compatibility=compatibility;
c.system=[c.K,sparse(compatibility);sparse(meanRow),sparse(1,1)];
c.factor=decomposition(c.system,'lu');
end
