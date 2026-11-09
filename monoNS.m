function [macro,info] = monoNS(macro,blks,cells,edges,gas,hot,nstep)
% Run monoNS for pure NS. GSIS calls the same solver with frozen HoT.
standalone=nargin==0;
if standalone
    fprintf('Starting pure NS: mesh and state initialization...\n'); drawnow;
    [blks,cells,edges,gas,macro]=nsw_init;
    hot.cell=zeros(blks.meshInfo(3),5); hot.boundary=zeros(numel(blks.far),4); nstep=4000;
end
N=blks.meshInfo(3); q=macro.qf; U=ns_conservative(q);
history=zeros(nstep,2);
if standalone; monitor=figure('Name','NSW pure NS monitor','NumberTitle','off'); end
timer=tic;
for it=1:nstep
    old=q;
    R=ns_residual(U,blks,cells,edges,gas,hot);
    dt=ns_timestep(q,macro.cfl,blks,cells,edges,gas);
    U=ns_implicit(U,R,dt,blks,cells,edges,gas,hot);
    q=ns_primitive(U);
    err=max(sqrt(sum((q-old).^2,1)./max(sum(old.^2,1),N))); history(it,:)=[err,max(abs(R),[],'all')];
    if standalone && (it==1 || mod(it,100)==0)
        fprintf('NS %d update %.6e equation %.6e elapsed %.1fs\n',it,err,history(it,2),toc(timer));
        ns_postplot(monitor,cells.xy(:,1),q,history(1:it,1));
    end
    if standalone && err<1e-8; break; end
end
macro.qf=q; macro.sigmaq=nsw_closure(q,blks,cells,edges,gas);
info.history=history(1:it,:); info.iterations=it; info.wall_time=toc(timer); final=ns_residual(U,blks,cells,edges,gas,hot); info.equation=max(abs(final),[],'all');
if standalone
    ns_postplot(monitor,cells.xy(:,1),q,info.history(:,1));
    out=fullfile(fileparts(mfilename('fullpath')),'results','unstructured'); if ~exist(out,'dir'); mkdir(out); end
    save(fullfile(out,'ns.mat'),'macro','blks','cells','edges','gas','info');
    writematrix([cells.xy,q,macro.sigmaq],fullfile(out,'ns_profile.dat'),'Delimiter','tab');
    writematrix([(1:it)',info.history],fullfile(out,'ns_history.dat'),'Delimiter','tab');
    fprintf('NS finished: %d, error %.8e, equation %.8e, %.2fs\n',it,err,info.equation,info.wall_time);
end
end

function ns_postplot(window,x,q,error)
if ~isgraphics(window); return; end
figure(window); clf(window); tiledlayout(2,2); [x,order]=sort(x);
columns=[1,2,4]; names={'rho','u_x','T'};
for j=1:3
    nexttile; plot(x,q(order,columns(j)),'LineWidth',1.2); xlabel('x'); ylabel(names{j}); xlim([-15,15]);
end
nexttile; loglog(1:numel(error),error,'LineWidth',1.2); xlabel('iteration'); ylabel('relative update'); grid on;
drawnow;
end

function U=ns_conservative(q)
U=[q(:,1),q(:,1).*q(:,2:3),q(:,1).*(1.5*q(:,4)+0.5*sum(q(:,2:3).^2,2))];
end

function q=ns_primitive(U)
q=[U(:,1),U(:,2:3)./U(:,1),(U(:,4)-0.5*sum(U(:,2:3).^2,2)./U(:,1))./(1.5*U(:,1))];
end

function dt=ns_timestep(q,cfl,blks,cells,edges,gas)
wave=abs(sum(q(blks.faces(:,3),2:3).*edges.enormal,2))+sqrt(gas.gamma*q(blks.faces(:,3),4));
dt=cfl*cells.A./(abs(blks.scatter)*(wave.*edges.elength));
end

function U=ns_implicit(U,R,dt,blks,cells,edges,gas,hot)
N=blks.meshInfo(3); row=[]; col=[]; value=[];
% Graph colouring groups columns with disjoint residual support.
for component=1:4
    for colour=1:max(cells.colour)
        group=find(cells.colour==colour); ids=group+(component-1)*N; epsU=1e-7*max(abs(U(ids)),1);
        Up=U; Up(ids)=Up(ids)+epsU; dR=ns_residual(Up,blks,cells,edges,gas,hot)-R;
        for j=1:numel(group)
            rr=reshape(cells.neighbours{group(j)}+(0:3)*N,[],1);
            row=[row;rr]; col=[col;repmat(ids(j),numel(rr),1)]; value=[value;dR(rr)/epsU(j)];
        end
    end
end
J=sparse(row,col,value,4*N,4*N); dU=(spdiags(repmat(1./dt,4,1),0,4*N,4*N)-J)\R(:);
U=U+reshape(dU,N,4);
end

function R=ns_residual(U,blks,cells,edges,gas,hot)
q=ns_primitive(U);
F=ns_convection(U,q,blks,edges,gas);
F=ns_viscous(F,q,blks,edges,gas,hot);
F=ns_boundary(F,q,blks,cells,edges,gas,hot);
R=-(blks.scatter*(F.*edges.elength))./cells.A;
end

function F=ns_convection(U,q,blks,edges,gas)
ids=blks.interior; L=blks.faces(ids,3); H=blks.faces(ids,4); n=edges.enormal(ids,:); a=q(L,:); b=q(H,:);
unL=sum(a(:,2:3).*n,2); unR=sum(b(:,2:3).*n,2); pL=a(:,1).*a(:,4); pR=b(:,1).*b(:,4);
FL=[U(L,1).*unL,U(L,2:3).*unL+pL.*n,(U(L,4)+pL).*unL];
FR=[U(H,1).*unR,U(H,2:3).*unR+pR.*n,(U(H,4)+pR).*unR];
wave=max(abs(unL)+sqrt(gas.gamma*a(:,4)),abs(unR)+sqrt(gas.gamma*b(:,4)));
F=zeros(blks.meshInfo(2),4); F(ids,:)=0.5*(FL+FR-wave.*(U(H,:)-U(L,:)));
end

function F=ns_viscous(F,q,blks,edges,gas,hot)
ids=blks.interior; L=blks.faces(ids,3); H=blks.faces(ids,4); n=edges.enormal(ids,:); a=q(L,:); b=q(H,:);
face=(a+b)/2; d=(b-a)./edges.distance(ids); mu=face(:,4).^gas.omega/gas.delta;
div=sum(d(:,2:3).*n,2); h=(hot.cell(L,:)+hot.cell(H,:))/2;
sxx=-mu.*(2*d(:,2).*n(:,1)-2/3*div)+h(:,1); sxy=-mu.*(d(:,2).*n(:,2)+d(:,3).*n(:,1))+h(:,2);
syy=-mu.*(2*d(:,3).*n(:,2)-2/3*div)+h(:,3); stress=[sxx.*n(:,1)+sxy.*n(:,2),sxy.*n(:,1)+syy.*n(:,2)];
heat=-2.5/gas.Pr*mu.*d(:,4)+sum(h(:,4:5).*n,2);
F(ids,:)=F(ids,:)+[zeros(numel(ids),1),stress,sum(stress.*face(:,2:3),2)+heat];
end

function F=ns_boundary(F,q,blks,cells,edges,gas,hot)
[~,boundary]=nsw_closure(q,blks,cells,edges,gas); F(blks.far,:)=boundary+hot.boundary;
ids=blks.sym; L=blks.faces(ids,3); n=edges.enormal(ids,:); F(ids,2:3)=q(L,1).*q(L,4).*n;
end
