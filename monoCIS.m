% First-order unstructured FV CIS; set demo_method='GSIS' to couple monoNS.
clearvars -except demo_method; clc;
if ~exist('demo_method','var'); demo_method='CIS'; end
is_gsis=strcmpi(demo_method,'GSIS'); clear demo_method;
if is_gsis; tag='GSIS'; else; tag='CIS'; end
fprintf('Starting %s: mesh and distribution initialization...\n',tag); drawnow;
[blks,cells,edges,gas,macro]=nsw_init;
N=blks.meshInfo(3); maxIter=6000; tolerance=1e-8; ns_inner_steps=100;
if is_gsis; tolerance=1e-6; end
[vx,vy]=meshgrid(linspace(-18,18,96),linspace(-15,15,80)); vx=vx(:)'; vy=vy(:)'; dv=36/95*30/79;
reflect=reshape(1:numel(vx),80,96); reflect=flipud(reflect); reflect=reflect(:)';
qf=macro.qf; f=qf(:,1)./(2*pi*qf(:,4)).*exp(-((vx-qf(:,2)).^2+(vy-qf(:,3)).^2)./(2*qf(:,4))); h=qf(:,4).*f;
qb=repmat(gas.left,numel(blks.boundary),1); right=blks.faces_tag(blks.boundary)==3; qb(right,:)=repmat(gas.right,nnz(right),1);
fBoundary=qb(:,1)./(2*pi*qb(:,4)).*exp(-((vx-qb(:,2)).^2+(vy-qb(:,3)).^2)./(2*qb(:,4))); hBoundary=qb(:,4).*fBoundary;
L=blks.faces(:,3); R=blks.faces_ghost(:,4); vn=edges.enormal(:,1)*vx+edges.enormal(:,2)*vy;
speed=vn.*edges.elength; pos=max(speed,0); neg=min(speed,0);
% Assemble outgoing coefficients directly from each cell's oriented faces.
outflow=zeros(N,numel(vx));
for i=1:N
    faces=cells.cell2face{i}; sign=2*(L(faces)==i)-1;
    outflow(i,:)=sum(max(sign.*speed(faces,:),0),1);
end
[qf,sigmaq]=dvm_moments(f,h,vx,vy,dv); total=sigmaq; active=false; history=zeros(maxIter,4);
monitor=figure('Name',['NSW ',tag,' monitor'],'NumberTitle','off'); timer=tic;
for iter=1:maxIter
    old=qf; nsold=macro.qf;
    if active; collision=macro.qf; heat=total(:,4:5); else; collision=qf; heat=sigmaq(:,4:5); end
    [G,H,nu]=dvm_collision(collision,heat,vx,vy,gas);
    [fext,hext]=dvm_boundary(f,h,fBoundary,hBoundary,blks,reflect);
    [fluxf,fluxh]=dvm_convection(fext,hext,L,R,pos,neg);
    rf=cells.A.*nu.*(G-f)-blks.scatter*fluxf; rh=cells.A.*nu.*(H-h)-blks.scatter*fluxh;
    [f,h]=dvm_implicit(f,h,rf,rh,nu,outflow,speed,blks,cells);
    [qf,sigmaq]=dvm_moments(f,h,vx,vy,dv);
    kerr=max(sqrt(sum((qf-old).^2,1)./max(sum(old.^2,1),N))); merr=0;
    if is_gsis && iter>=10
        hot=dvm_hot(f,h,qf,sigmaq,fBoundary,hBoundary,vn,vx,vy,dv,blks,cells,edges,gas);
        if iter==10
            macro.qf=qf; total=sigmaq; active=true;
        else
            [macro,ns_info]=monoNS(macro,blks,cells,edges,gas,hot,ns_inner_steps);
            total=macro.sigmaq+hot.cell; merr=max(sqrt(sum((macro.qf-nsold).^2,1)./max(sum(nsold.^2,1),N)));
        end
    end
    err=max(kerr,merr); history(iter,:)=[kerr,merr,err,toc(timer)];
    if iter==1 || mod(iter,5)==0 || (is_gsis && mod(iter,10)==0)
        fprintf('%s %d kinetic %.6e macro %.6e elapsed %.1fs\n',tag,iter,kerr,merr,history(iter,4)); drawnow limitrate;
    end
    if iter==1 || mod(iter,100)==0 || (is_gsis && mod(iter,10)==0); dvm_postplot(monitor,cells.xy(:,1),qf,history(1:iter,3)); end
    if err<tolerance && (~is_gsis || iter>10); break; end
end
history=history(1:iter,:); out=fullfile(fileparts(mfilename('fullpath')),'results','unstructured'); if ~exist(out,'dir'); mkdir(out); end
dvm_postplot(monitor,cells.xy(:,1),qf,history(:,3));
if is_gsis; tag='gsis'; else; tag='cis'; end
save(fullfile(out,[tag,'.mat']),'f','h','qf','sigmaq','macro','blks','cells','edges','gas','vx','vy','dv','history','iter');
if is_gsis; save(fullfile(out,[tag,'.mat']),'hot','total','ns_inner_steps','ns_info','-append'); end
writematrix([cells.xy,qf,sigmaq],fullfile(out,[tag,'_profile.dat']),'Delimiter','tab');
writematrix([(1:iter)',history],fullfile(out,[tag,'_history.dat']),'Delimiter','tab');
fprintf('%s finished: %d error %.8e time %.2fs\n',tag,iter,err,history(end,4));

function [G,H,nu]=dvm_collision(q,heat,vx,vy,gas)
rho=q(:,1); T=q(:,4); cx=vx-q(:,2); cy=vy-q(:,3); c2=cx.^2+cy.^2;
E=rho./(2*pi*T).*exp(-c2./(2*T)); A=(1-gas.Pr)./(5*rho.*T.^2).*(heat(:,1).*cx+heat(:,2).*cy);
G=E.*(1+A.*(c2./T-4)); H=T.*E.*(1+A.*(c2./T-2)); nu=gas.delta*rho.*T.^(1-gas.omega);
end

function dvm_postplot(window,x,q,error)
if ~isgraphics(window); return; end
figure(window); clf(window); tiledlayout(2,2); [x,order]=sort(x);
columns=[1,2,4]; names={'rho','u_x','T'};
for j=1:3
    nexttile; plot(x,q(order,columns(j)),'LineWidth',1.2); xlabel('x'); ylabel(names{j}); xlim([-15,15]);
end
nexttile; loglog(1:numel(error),error,'LineWidth',1.2); xlabel('iteration'); ylabel('relative update'); grid on;
drawnow;
end

function [fext,hext]=dvm_boundary(f,h,fb,hb,blks,reflect)
L=blks.faces(:,3); sym=blks.faces_tag(blks.boundary)==1;
fb(sym,:)=f(L(blks.sym),reflect); hb(sym,:)=h(L(blks.sym),reflect);
fext=[f;fb]; hext=[h;hb];
end

function [fluxf,fluxh]=dvm_convection(fext,hext,L,R,pos,neg)
fluxf=pos.*fext(L,:)+neg.*fext(R,:); fluxh=pos.*hext(L,:)+neg.*hext(R,:);
end

function [f,h]=dvm_implicit(f,h,rf,rh,nu,outflow,speed,blks,cells)
N=blks.meshInfo(3); L=blks.faces(:,3);
D=cells.A.*(nu+1e-5)+outflow; df=zeros(size(f)); dh=zeros(size(h));
for pass=1:2
    if pass==1; sequence=1:N; else; sequence=N:-1:1; end
    for i=sequence
        a=rf(i,:); b=rh(i,:);
        for face=cells.cell2face{i}
            if L(face)==i; j=blks.faces(face,4); s=speed(face,:); else; j=L(face); s=-speed(face,:); end
            if j>0; a=a-min(s,0).*df(j,:); b=b-min(s,0).*dh(j,:); end
        end
        df(i,:)=a./D(i,:); dh(i,:)=b./D(i,:);
    end
end
f=f+df; h=h+dh;
end

function hot=dvm_hot(f,h,q,sigmaq,fb,hb,vn,vx,vy,dv,blks,cells,edges,gas)
[closure,Fns]=nsw_closure(q,blks,cells,edges,gas); hot.cell=sigmaq-closure;
L=blks.faces(:,3); R=blks.faces_ghost(:,4); fext=[f;fb]; hext=[h;hb]; ids=blks.far;
ff=max(vn(ids,:),0).*fext(L(ids),:)+min(vn(ids,:),0).*fext(R(ids),:);
hh=max(vn(ids,:),0).*hext(L(ids),:)+min(vn(ids,:),0).*hext(R(ids),:);
Fkin=dv*[sum(ff,2),sum(ff.*vx,2),sum(ff.*vy,2),0.5*sum((vx.^2+vy.^2).*ff+hh,2)];
hot.boundary=Fkin-Fns;
end

function [q,sq]=dvm_moments(f,h,vx,vy,dv)
rho=sum(f,2)*dv; u=sum(f.*vx,2)*dv./rho; v=sum(f.*vy,2)*dv./rho; cx=vx-u; cy=vy-v;
T=sum(f.*(cx.^2+cy.^2)+h,2)*dv./(3*rho); energy=(cx.^2+cy.^2).*f+h;
q=[rho,u,v,T]; sq=dv*[sum(f.*cx.^2,2),sum(f.*cx.*cy,2),sum(f.*cy.^2,2),0.5*sum(cx.*energy,2),0.5*sum(cy.*energy,2)];
sq(:,[1,3])=sq(:,[1,3])-rho.*T;
end
