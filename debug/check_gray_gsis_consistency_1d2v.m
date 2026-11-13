function check_gray_gsis_consistency_1d2v
% Independent wall tangent, seven kinetic weak moments and monolithic solve.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
old=warning('off','Channel1D2V:NotConverged');restore=onCleanup(@()warning(old));
k=mono1D2V(struct('flow','Poiseuille','cells',4,'delta',3, ...
    'normalQuadrature','hermite','velocityOrder',[8 8],'maxIter',1, ...
    'plot',false,'save',false,'reference',false));
k.mesh.edges=[-.5;-.37;-.05;.21;.5];
v=k.velocity;N=4;nv=numel(v.weight);rng(20261110);
k.H=.01*randn(2*N,nv); % Jumps and all kinetic moments, not just shear parity.
meanH=sum(diff(k.mesh.edges).*(k.H(1:2:end,:)+k.H(2:2:end,:))*v.weight/2);
k.H=k.H-meanH/sum(v.weight);
[~,~,h,hot]=gsis_macro_from_kinetic_1d2v(k);
bc=hot.boundary;c=h.cache;

% Finite difference actual discrete diffuse reflection, including n_wall.
D=[1 0 0 0;0 1 0 0;0 0 1 0;-1 0 0 1];
step=1e-5;
for component=1:4
    change=D(component,:)/v.gram*v.basis';
    plus=k;minus=k;plus.H=k.H+step*change;minus.H=k.H-step*change;
    bp=kinetic_macro_wall_bc_1d2v(plus);bm=kinetic_macro_wall_bc_1d2v(minus);
    fd=(bp.flux-bm.flux)/(2*step);
    exact=squeeze(bc.jacobian(:,component,:))';
    assert(max(abs(fd-exact),[],'all')<1e-10);
end
assert(max(abs(bc.flux(:,1)))<1e-14);
assert(max(abs(bc.jacobian(1,:,:)),[],'all')<1e-14);
fprintf('PASS all four wall Jacobian columns vs actual diffuse-flux finite differences.\n');

% Form the full Nv-dimensional DG residual first, then integrate its moments.
test=[ones(nv,1),2*v.vx,2*v.vy,v.speed2-2, ...
    v.vx.^2-v.vy.^2,2*v.vx.*v.vy,v.vx.*(v.speed2-2)];
q=k.H*v.moment;gain=(q/v.gram)*v.basis';
for cellId=1:N
    nodes=2*cellId-1:2*cellId;
    R=(-c.G*k.H(nodes,:)+[-hot.kineticTrace(cellId,:); ...
        hot.kineticTrace(cellId+1,:)]).*v.vx' ...
        +c.M(:,:,cellId)*(k.config.delta*(k.H(nodes,:)-gain(nodes,:)) ...
        -2*k.config.acceleration*v.vy');
    moments=R*(v.weight.*test);
    assert(max(abs(moments-hot.consistency.kineticWeakMomentResidual(nodes,:)),[],'all')<1e-13);
    z=hot.consistency.baseFields(:,cellId);
    trace=reshape(hot.consistency.baseTrace(cellId:cellId+1,:)',[],1);
    macroResidual=c.A(:,:,cellId)*z-c.B(:,:,cellId)*trace-h.localRhs(:,cellId);
    target=reshape([moments(:,1:4),-moments(:,5:7)],[],1);
    assert(max(abs(macroResidual-target))<1e-13);
end
fprintf('PASS seven macro weak residuals vs moments of the full kinetic DG residual.\n');

% Independently assemble uncondensed fields + traces + density multiplier.
nf=14*N;nt=4*(N+1);A=zeros(nf+nt+1);b=zeros(nf+nt+1,1);
for i=1:N
    fields=14*(i-1)+(1:14);faces=nf+c.traceIds(:,i);
    A(fields,fields)=c.A(:,:,i);A(fields,faces)=-c.B(:,:,i);
    A(faces,fields)=A(faces,fields)+c.C(:,:,i);
    A(faces,faces)=A(faces,faces)+c.D(:,:,i);b(fields)=h.localRhs(:,i);
    measure=sum(c.M(:,:,i),1)/(c.edges(end)-c.edges(1));
    A(end,fields(1:2))=measure;A(end,fields(7:8))=-measure;
end
for side=1:2
    i=1+(side-1)*(N-1);localRows=(side-1)*4+(1:4);
    rows=nf+c.traceIds(localRows,i);fields=14*(i-1)+(1:14);
    E=zeros(4,14);for j=1:4;E(j,2*(j-1)+side)=1;end
    A(rows,fields)=A(rows,fields)-bc.jacobian(:,:,side)*E;
    b(rows)=bc.flux(side,:)'-bc.jacobian(:,:,side)*bc.state(side,:)';
end
A(nf+(1:4:nt),end)=1/sqrt(N+1);b(end)=bc.meanDensity;
z=A\b;
assert(max(abs(A*z-b))<1e-12);
assert(max(abs(reshape(z(1:nf),14,N)-h.localFields),[],'all')<1e-11);
assert(max(abs(reshape(z(nf+(1:nt)),4,[])'-h.Uhat),[],'all')<1e-11);
assert(abs(z(end)-h.bc.gaugeMultiplier)<1e-12);
assert(max([h.boundaryFluxResidual,h.bc.densityConstraintResidual, ...
    h.traceResidual,h.localResidual,h.fluxContinuityResidual,h.conservationResidual])<1e-12);
fprintf('PASS affine-wall/density-gauge condensed solve vs independent monolithic system.\n');
fprintf('All Gray-adapted BGK consistency checks passed.\n');
end
