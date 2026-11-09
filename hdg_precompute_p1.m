function cache = hdg_precompute_p1(mesh,delta,param)
%HDG_PRECOMPUTE_P1 Local P1 mixed BGK operators and condensed trace matrix.
% mesh.edges: strictly increasing cell faces (nonuniform meshes allowed).
% Reference-code pressure form, adapted from 2D3V Shakhov to 1D2V BGK.
% param.stabilization: positive [s_p s_ux s_uy s_T], default [1 delta delta delta].
% Local ordering: [p(1:2),ux(1:2),uy(1:2),T(1:2),Lx(1:2),Ly(1:2),E(1:2)].
% p=rho+T. Trace ordering: [pL uxL uyL TL pR uxR uyR TR].
% A*U=f+B*Uhat; outward numerical fluxes = C*U+D*Uhat.
if nargin<3; param=struct; end
assert(isstruct(mesh) && isfield(mesh,'edges'),'mesh.edges is required.');
edges=mesh.edges(:);
validateattributes(edges,{'numeric'},{'vector','real','finite','increasing'});
assert(numel(edges)>=2,'At least one cell is required.');
validateattributes(delta,{'numeric'},{'scalar','real','finite','positive'});
s=[1 delta delta delta];
if isfield(param,'stabilization'); s=param.stabilization; end
validateattributes(s,{'numeric'},{'vector','numel',4,'real','finite','positive'});
s=s(:)'; N=numel(edges)-1; lengthCell=diff(edges);
cache=struct('edges',edges,'count',N,'delta',delta,'stabilization',s, ...
    'formulation','pressure-reference-v1', ...
    'nodes',reshape([edges(1:end-1),edges(2:end)]',[],1), ...
    'G',[-1 -1;1 1]/2,'M',zeros(2,2,N),'A',zeros(14,14,N), ...
    'B',zeros(14,8,N),'C',zeros(8,14,N),'D',zeros(8,8,N), ...
    'response',zeros(14,8,N),'traceIds',zeros(8,N),'localRcond',zeros(N,1));
% Scaled momentum rows, and heat equation after subtracting mass balance:
% F=[ux; p-Lx/delta; -Ly/delta; -E/delta].
F=zeros(4,7); F(1,2)=1; F(2,[1 5])=[1 -1/delta];
F(3,6)=-1/delta; F(4,7)=-1/delta;
% Mass uses LOCAL ux; pressure in momentum uses the common trace.
localFace=F; localFace(2,1)=0;
traceFace=zeros(4); traceFace(2,1)=1;
select=zeros(3,7); select(:,2:4)=eye(3);
rows=zeros(64,N); cols=rows; values=rows;
for i=1:N
    M=lengthCell(i)*[2 1;1 2]/6; G=cache.G;
    A=[-kron(F,G);kron(select,G)];
    A(9:14,9:14)=A(9:14,9:14)+kron(eye(3),M);
    B=zeros(14,8); C=zeros(8,14); D=zeros(8,8);
    for side=1:2
        b=zeros(2,1); b(side)=1; normal=2*side-3;
        faceRows=(side-1)*4+(1:4);
        Clocal=kron(normal*localFace+[diag(s),zeros(4,3)],b');
        Dlocal=normal*traceFace-diag(s);
        A(1:8,:)=A(1:8,:)+kron(eye(4),b)*Clocal;
        B(1:8,faceRows)=-kron(eye(4),b)*Dlocal;
        B(9:14,faceRows)=normal*kron(select(:,1:4),b);
        C(faceRows,:)=Clocal; D(faceRows,faceRows)=Dlocal;
    end
    cache.localRcond(i)=rcond(A);
    assert(cache.localRcond(i)>eps,'Singular local HDG matrix in cell %d.',i);
    response=A\B; Schur=C*response+D;
    ids=4*(i-1)+(1:8)'; [rr,cc]=ndgrid(ids,ids);
    rows(:,i)=rr(:); cols(:,i)=cc(:); values(:,i)=Schur(:);
    cache.M(:,:,i)=M; cache.A(:,:,i)=A; cache.B(:,:,i)=B;
    cache.C(:,:,i)=C; cache.D(:,:,i)=D;
    cache.response(:,:,i)=response; cache.traceIds(:,i)=ids;
end
cache.K=sparse(rows(:),cols(:),values(:),4*(N+1),4*(N+1));
end
