function check_gsis_cavity_2d2v
% Independent quadrature weak form, monolithic recovery, and no-lift checks.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
old=warning('off','Cavity:NotConverged');cleanup=onCleanup(@()warning(old));
cfg=struct('cells',[3 2],'meshStretch',1.3,'velocityOrder',[8 8], ...
    'quadrature','hermite','delta',3,'maxIter',1,'plot',false,'save',false, ...
    'macroParam',struct('kineticConsistency',false));
k=monoCavity2D2V(cfg);v=k.velocity;N=k.mesh.count;
c=hdg_precompute_q1_2d2v(k.mesh,v,cfg.delta);
[macro,h]=gsis_macro_from_kinetic_2d2v(k,struct('cache',c));
nf=32*N;nt=8*c.faceCount;A=zeros(nf+nt+1);b=zeros(nf+nt+1,1);error=0;
for cellId=1:N
    dx=k.mesh.dx(cellId);dy=k.mesh.dy(cellId);
    Ar=zeros(32);Br=Ar;Cr=Ar;Dr=Ar;
    for component=1:32
        z=zeros(32,1);z(component)=1;
        [r,f]=weak(z,zeros(32,1),dx,dy,cfg.delta,c.stabilization);
        Ar(:,component)=r;Cr(:,component)=f;
        [r,f]=weak(zeros(32,1),z,dx,dy,cfg.delta,c.stabilization);
        Br(:,component)=-r;Dr(:,component)=f;
    end
    error=max(error,max(abs([Ar-c.A(:,:,cellId),Br-c.B(:,:,cellId)]),[],'all'));
    rows=32*(cellId-1)+(1:32);traces=nf+c.ids(:,cellId);
    A(rows,rows)=Ar;A(rows,traces)=-Br;b(rows)=h.localRhs(:,cellId);
    A(traces,rows)=A(traces,rows)+Cr;A(traces,traces)=A(traces,traces)+Dr;
    M=dx*dy*[4 2 2 1;2 4 1 2;2 1 4 2;1 2 2 4]/36;
    A(end,rows(1:4))=sum(M,1);A(end,rows(13:16))=-sum(M,1);
    for side=1:4
        localRows=8*(side-1)+(1:8);
        error=max(error,max(abs(Cr(localRows,:)-c.C(:,:,side,cellId)),[],'all'));
        error=max(error,max(abs(Dr(localRows,localRows)-c.D(:,:,side,cellId)),[],'all'));
        if k.mesh.neighbor(cellId,side)>0;continue;end
        E=c.E{side};length=dy;if side>2;length=dx;end
        Mf=length*[2 1;1 2]/6;J=c.wallJ(:,:,side);
        A(traces(localRows),rows(1:16))=A(traces(localRows),rows(1:16))-kron(J,Mf*E');
        [flux,state]=wall_flux(k,side,cellId);
        target=Mf*(flux-state*J');b(traces(localRows))=target(:);
        % Finite difference the actual diffuse map for all primitive columns.
        D=[1 0 0 0;0 1 0 0;0 0 1 0;-1 0 0 1];
        for component=1:4
            deltaH=D(component,:)/v.gram*v.basis';epsH=1e-5;
            kp=k;km=k;kp.H=k.H+epsH*deltaH;km.H=k.H-epsH*deltaH;
            fp=wall_flux(kp,side,cellId);fm=wall_flux(km,side,cellId);
            assert(max(abs((fp-fm)/(2*epsH)-J(:,component)'),[],'all')<1e-10);
        end
    end
    % Explicitly verify component/node/cell ordering of macro feedback.
    for component=1:4
        nodes=4*(cellId-1)+(1:4);
        assert(max(abs(macro.pressureNodal(nodes,component)-h.localFields(4*(component-1)+(1:4),cellId)))<1e-14);
    end
end
assert(error<1e-12);
A(nf+(1:nt),end)=c.compatibility;sol=A\b;
assert(max(abs(A*sol-b))<1e-11);
assert(max(abs(reshape(sol(1:nf),32,N)-h.localFields),[],'all')<1e-10);
assert(max(abs(sol(nf+(1:nt))-h.trace))<1e-10);
assert(max([h.localResidual,h.traceResidual,h.wallFluxResidual,abs(h.densityGaugeMultiplier)])<1e-11);
fprintf('PASS independent 2D Q1 weak integration, all face normals, diffuse Jacobians, condensed vs monolithic recovery; operator error %.3e.\n',error);
g=monoGSISCavity2D2V(cfg);
assert(isequal(k.H,g.H),'Macro solve must never lift or overwrite H.');
assert(max(abs(g.synthetic.sourceMoments-g.synthetic.macro.nodal),[],'all')<1e-12);
assert(~g.synthetic.hdg.internalResidualCorrection);
fprintf('PASS identical first kinetic step, full macro source return, no H lift or internal residual load correction.\n');

% Optional consistency load: independently form the full Nv-dimensional
% kinetic residual first, then integrate its eight velocity moments.
[~,hc]=gsis_macro_from_kinetic_2d2v(k,struct('cache',c,'kineticConsistency',true));
test=[ones(size(v.vx)),2*v.vx,2*v.vy,v.speed2-2,v.vx.^2-v.vy.^2, ...
    2*v.vx.*v.vy,v.vx.*(v.speed2-2),v.vy.*(v.speed2-2)];
gain=((k.H*v.moment)/v.gram)*v.basis';nodeIds={[1 3],[2 4],[1 2],[3 4]};opposite=[2 1 4 3];
bc=b;
for cellId=1:N
    ids=4*(cellId-1)+(1:4);R=(-c.Gx(:,:,cellId)*k.H(ids,:)).*v.vx' ...
        +(-c.Gy(:,:,cellId)*k.H(ids,:)).*v.vy' ...
        +cfg.delta*c.M(:,:,cellId)*(k.H(ids,:)-gain(ids,:));
    for side=1:4
        normal=c.normals(side,:);cn=normal(1)*v.vx+normal(2)*v.vy;
        nodes=4*(cellId-1)+nodeIds{side};trace=k.H(nodes,:);in=cn<0;
        if k.mesh.neighbor(cellId,side)>0
            other=4*(k.mesh.neighbor(cellId,side)-1)+nodeIds{opposite(side)};
            trace(:,in)=k.H(other,in);
        else
            [~,~,trace]=wall_flux(k,side,cellId);
        end
        R=R+(c.E{side}*c.faceM(:,:,side,cellId)*trace).*cn';
    end
    moments=R*(v.weight.*test);target=reshape([moments(:,1:4),-moments(:,5:8)],[],1);
    assert(max(abs(target-hc.targetKineticMomentResidual(:,cellId)))<1e-13);
    hat0=reshape(hc.baseTrace(:,c.faces(cellId,:)),[],1);
    actual=c.A(:,:,cellId)*hc.baseFields(:,cellId)-c.B(:,:,cellId)*hat0-hc.localRhs(:,cellId);
    assert(max(abs(actual-target))<1e-13);
    bc(32*(cellId-1)+(1:32))=hc.localRhs(:,cellId);
end
corrected=A\bc;
assert(max(abs(reshape(corrected(1:nf),32,N)-hc.localFields),[],'all')<1e-10);
assert(max(abs(corrected(nf+(1:nt))-hc.trace))<1e-10);
cfg=struct('cells',4,'velocityOrder',[8 8],'quadrature','hermite','delta',10, ...
    'tol',1e-10,'plot',false,'save',false,'macroParam',struct('kineticConsistency',true));
cis=monoCavity2D2V(cfg);[m,~]=gsis_macro_from_kinetic_2d2v(cis,cfg.macroParam);
assert(norm(m.nodal-cis.nodalMoments,'fro')/norm(cis.nodalMoments,'fro')<1e-7);
cfg.tol=1e-8;g=monoGSISCavity2D2V(cfg);
assert(g.converged && g.diagnostics.originalConverged);
assert(norm(g.nodalMoments-cis.nodalMoments,'fro')/norm(cis.nodalMoments,'fro')<1e-6);
fprintf('PASS optional eight-moment consistency vs full kinetic residual, monolithic solve, one-step CIS fixed point and coupled CIS agreement.\n');
cfg.macroParam.relaxation=.1;g=monoGSISCavity2D2V(cfg);
expected=g.nodalMoments+.1*(g.synthetic.macro.nodal-g.nodalMoments);
measure=repelem(g.mesh.area/4,4);expected(:,1)=expected(:,1)-measure'*expected(:,1);
assert(max(abs(expected-g.synthetic.sourceMoments),[],'all')<1e-12);
assert(g.converged && g.diagnostics.originalConverged && g.synthetic.relaxation==.1);
assert(norm(g.nodalMoments-cis.nodalMoments,'fro')/norm(cis.nodalMoments,'fro')<1e-6);
fprintf('PASS relaxed macro feedback formula, original-equation convergence and unchanged CIS fixed point.\n');
end

function [flux,state,trace]=wall_flux(k,side,cellId)
v=k.velocity;nodes={[1 3],[2 4],[1 2],[3 4]};normals=[-1 0;1 0;0 -1;0 1];
ids=4*(cellId-1)+nodes{side};trace=k.H(ids,:);cn=normals(side,1)*v.vx+normals(side,2)*v.vy;
in=cn<0;out=cn>0;motion=zeros(1,numel(cn));if side==4;motion=2*k.config.lidSpeed*v.vx';end
rho=-(trace(:,out)*(v.weight(out).*cn(out))+sum(motion(in)'.*v.weight(in).*cn(in)))/sum(v.weight(in).*cn(in));
trace(:,in)=rho+motion(in);
flux=trace*((v.weight.*cn).*[ones(size(cn)),2*v.vx,2*v.vy,v.speed2-2]);
q=k.H(ids,:)*v.moment;state=[q(:,1)+q(:,4),q(:,2:4)];
end

function [residual,faceResidual]=weak(z,hat,dx,dy,delta,s)
Z=reshape(z,4,8);hat=reshape(hat,2,4,4);residual=zeros(4,8);faceResidual=zeros(2,4,4);
gp=(1+[-1 1]/sqrt(3))/2;
for x=gp
    for y=gp
        [p,px,py]=shape(x,y,dx,dy);u=p'*Z;
        Fx=[u(2),u(1)-u(5)/delta,-u(6)/delta,-u(7)/delta];
        Fy=[u(3),-u(6)/delta,u(1)+u(5)/delta,-u(8)/delta];
        residual(:,1:4)=residual(:,1:4)-dx*dy/4*(px*Fx+py*Fy);
        residual(:,5:8)=residual(:,5:8)+dx*dy/4*(p*u(5:8) ...
            +px*[u(2),u(3),u(4),0]+py*[-u(3),u(2),0,u(4)]);
    end
end
normals=[-1 0;1 0;0 -1;0 1];
for side=1:4
    n=normals(side,:);length=dy;if side>2;length=dx;end
    for t=gp
        if side==1;p=shape(0,t,dx,dy);elseif side==2;p=shape(1,t,dx,dy);
        elseif side==3;p=shape(t,0,dx,dy);else;p=shape(t,1,dx,dy);end
        edge=[1-t;t];u=p'*Z;h=edge'*hat(:,:,side);
        f=[n(1)*u(2)+n(2)*u(3),n(1)*(h(1)-u(5)/delta)-n(2)*u(6)/delta, ...
            -n(1)*u(6)/delta+n(2)*(h(1)+u(5)/delta),-(n(1)*u(7)+n(2)*u(8))/delta] ...
            +s.*(u(1:4)-h);
        residual(:,1:4)=residual(:,1:4)+length/2*p*f;
        residual(:,5:8)=residual(:,5:8)-length/2*p* ...
            [n(1)*h(2)-n(2)*h(3),n(2)*h(2)+n(1)*h(3),n(1)*h(4),n(2)*h(4)];
        faceResidual(:,:,side)=faceResidual(:,:,side)+length/2*edge*f;
    end
end
residual=residual(:);faceResidual=faceResidual(:);
end

function [p,px,py]=shape(x,y,dx,dy)
p=[(1-x)*(1-y);x*(1-y);(1-x)*y;x*y];
px=[-(1-y);1-y;-y;y]/dx;py=[-(1-x);-x;1-x;x]/dy;
end
