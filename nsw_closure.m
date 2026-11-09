function [sigmaq,boundaryFlux] = nsw_closure(q,blks,cells,edges,gas)
% Physical NS closure and outward half-Maxwell far-field flux.
L=blks.faces(:,3); R=blks.faces_ghost(:,4); b=blks.boundary;
ghost=q(L(b),:); tag=blks.faces_tag(b);
ghost(tag==2,:)=repmat(gas.left,nnz(tag==2),1); ghost(tag==3,:)=repmat(gas.right,nnz(tag==3),1);
qext=[q;ghost]; face=(qext(L,:)+qext(R,:))/2;
gx=(blks.scatter*(face.*edges.enormal(:,1).*edges.elength))./cells.A;
gy=(blks.scatter*(face.*edges.enormal(:,2).*edges.elength))./cells.A;
mu=q(:,4).^gas.omega/gas.delta; div=gx(:,2)+gy(:,3);
sigmaq=[-mu.*(2*gx(:,2)-2/3*div),-mu.*(gx(:,3)+gy(:,2)),-mu.*(2*gy(:,3)-2/3*div),-2.5/gas.Pr*mu.*gx(:,4),-2.5/gas.Pr*mu.*gy(:,4)];
ids=blks.far; n=edges.enormal(ids,:); a=q(L(ids),:); z=repmat(gas.left,numel(ids),1); z(blks.faces_tag(ids)==3,:)=repmat(gas.right,nnz(blks.faces_tag(ids)==3),1);
boundaryFlux=half(a,n,1)+half(z,n,-1);
end

function F=half(q,n,sgn)
rho=q(:,1); T=q(:,4); un=sum(q(:,2:3).*n,2); ut=-q(:,2).*n(:,2)+q(:,3).*n(:,1);
m0=0.5*erfc(-sgn*un./sqrt(2*T)); m1=un.*m0+sgn*sqrt(T/(2*pi)).*exp(-un.^2./(2*T));
m2=un.*m1+T.*m0; m3=un.*m2+2*T.*m1;
F=rho.*[m1,m2.*n(:,1)-ut.*m1.*n(:,2),m2.*n(:,2)+ut.*m1.*n(:,1),0.5*(m3+(ut.^2+2*T).*m1)];
end
