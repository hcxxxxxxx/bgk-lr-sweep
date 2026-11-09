function [blks,cells,edges,gas,macro] = nsw_init
% Fluent topology is retained; coordinates are scaled to the demo domain.
root=fileparts(mfilename('fullpath'));
blks=importFluentMesh2(fullfile(root,'mesh','nsw400.cas'),1);
x=blks.points(:,1); y=blks.points(:,2);
blks.points(:,1)=120*(x-min(x))/(max(x)-min(x))-60;
blks.points(:,2)=(y-min(y))/(max(y)-min(y));
N=blks.meshInfo(3); F=blks.meshInfo(2); faces=blks.faces;
cells.xy=zeros(N,2); cells.A=zeros(N,1); cells.cell2face=cell(N,1);
for k=1:F
    cells.cell2face{faces(k,3)}(end+1)=k;
    if faces(k,4)>0; cells.cell2face{faces(k,4)}(end+1)=k; end
end
for i=1:N
    nodes=unique(reshape(faces(cells.cell2face{i},1:2),[],1)); xy=blks.points(nodes,:); center=mean(xy,1);
    [~,j]=sort(atan2(xy(:,2)-center(2),xy(:,1)-center(1))); xy=xy(j,:); xy1=xy([2:end,1],:);
    a=xy(:,1).*xy1(:,2)-xy1(:,1).*xy(:,2); cells.A(i)=sum(a)/2;
    cells.xy(i,:)=sum((xy+xy1).*a,1)/(6*cells.A(i));
end
v=blks.points(faces(:,2),:)-blks.points(faces(:,1),:);
edges.xy=(blks.points(faces(:,2),:)+blks.points(faces(:,1),:))/2;
edges.elength=hypot(v(:,1),v(:,2)); edges.enormal=[v(:,2),-v(:,1)]./edges.elength;
flip=sum(edges.enormal.*(edges.xy-cells.xy(faces(:,3),:)),2)<0; edges.enormal(flip,:)=-edges.enormal(flip,:);
b=find(faces(:,4)==0); blks.faces_ghost=faces; blks.faces_ghost(b,4)=N+(1:numel(b));
cells.xy_ghost=[cells.xy;2*edges.xy(b,:)-cells.xy(faces(b,3),:)];
edges.distance=sum((cells.xy_ghost(blks.faces_ghost(:,4),:)-cells.xy(faces(:,3),:)).*edges.enormal,2);
internal=find(faces(:,4)>0); blks.interior=internal; blks.boundary=b;
blks.faces_tag=zeros(F,1); blks.faces_tag(b)=1;
blks.faces_tag(b(edges.xy(b,1)<-59.999))=2; blks.faces_tag(b(edges.xy(b,1)>59.999))=3;
blks.far=find(blks.faces_tag>=2); blks.sym=find(blks.faces_tag==1);
blks.scatter=sparse([faces(:,3);faces(internal,4)],[(1:F)';internal],[ones(F,1);-ones(numel(internal),1)],N,F);
adj=spones(sparse(faces(internal,3),faces(internal,4),1,N,N)); adj=spones(adj+adj'+speye(N));
conflict=spones(adj*adj); cells.colour=zeros(N,1); cells.neighbours=cell(N,1);
for i=1:N
    used=cells.colour(find(conflict(i,:))); colour=1; while any(used==colour); colour=colour+1; end
    cells.colour(i)=colour; cells.neighbours{i}=find(adj(:,i));
end
gas=struct('gamma',5/3,'Pr',2/3,'omega',0.74,'delta',sqrt(pi/2),'Ma',3);
r=(gas.gamma+1)*gas.Ma^2/(2+(gas.gamma-1)*gas.Ma^2);
t=(2*gas.gamma*gas.Ma^2-gas.gamma+1)*((gas.gamma-1)*gas.Ma^2+2)/((gas.gamma+1)^2*gas.Ma^2);
gas.left=[1,gas.Ma*sqrt(gas.gamma),0,1]; gas.right=[r,gas.left(2)/r,0,t];
macro.qf=repmat(gas.left,N,1); right=cells.xy(:,1)>=0; macro.qf(right,:)=repmat(gas.right,nnz(right),1); macro.cfl=5;
assert(N>=50 && all(cells.A>0) && all(edges.distance>0));
assert(all(abs(edges.enormal(blks.sym,1))<1e-10),'Demo reflection requires horizontal symmetry faces');
assert(nnz(blks.faces_tag==2)>0 && nnz(blks.faces_tag==3)>0,'Missing upstream/downstream boundary');
assert(max(abs(blks.scatter*(edges.enormal.*edges.elength)),[],'all')<1e-10);
fprintf('Fluent FV: %d cells, %d faces, x=[-60,60], Ma=3\n',N,F);
end
