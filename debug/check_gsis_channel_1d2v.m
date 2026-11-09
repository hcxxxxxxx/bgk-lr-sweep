function check_gsis_channel_1d2v
%CHECK_GSIS_CHANNEL_1D2V Source-only coupling and reference-style data checks.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
diary(fullfile(root,'debug','gsis_channel_validation.log'));
cleanup=onCleanup(@() diary('off'));
base=struct('flow','Poiseuille','delta',10,'cells',8,'velocityOrder',[32 8], ...
    'plot',false,'save',false,'reference',false,'maxIter',1,'reportEvery',1000, ...
    'stopCriterion','equation','tol',1e-10);
old=warning('off','Channel1D2V:NotConverged'); restore=onCleanup(@() warning(old));
c=mono1D2V(base); g=monoGSIS1D2V(base);
% First kinetic sweep must be identical: the macro solve never lifts H.
assert(isequal(c.H,g.H) && isequal(c.kineticTrace,g.kineticTrace));
assert(isequal(g.synthetic.hdg.cache.stabilization,[1 10 10 10]));
expected=g.nodalMoments+g.synthetic.omega.*(g.synthetic.macro.nodal-g.nodalMoments);
expected(:,1)=expected(:,1)-mean(expected(:,1));
assert(max(abs(expected-g.synthetic.sourceMoments),[],'all')<1e-14);
assert(norm(g.synthetic.sourceMoments-g.nodalMoments,'fro')>1e-6);
fprintf('PASS identical first kinetic sweep, no H lift, mixed source and unchanged stabilization.\n');

% Independently solve the SECOND kinetic step using the saved mixed source.
two=base; two.maxIter=2; r=monoGSIS1D2V(two);
v=g.velocity; N=g.mesh.count; M=g.mesh.dx*[2 1;1 2]/6; G=[-1 -1;1 1]/2;
source=(g.synthetic.sourceMoments/v.gram)*v.basis';
rhs=kron(speye(N),M)*(base.delta*source+2*g.config.acceleration*v.vy');
wall=zeros(2,numel(v.vx));
wall(1,:)=g.H(1,v.vx<0)*(v.weight(v.vx<0).*abs(v.vx(v.vx<0)))/v.halfFlux(1);
wall(2,:)=g.H(end,v.vx>0)*(v.weight(v.vx>0).*abs(v.vx(v.vx>0)))/v.halfFlux(2);
wall=wall+2*g.config.wallVelocity(:)*v.vy';
rhs(1,v.inLeft)=rhs(1,v.inLeft)+wall(1,v.inLeft).*v.vx(v.inLeft)';
rhs(end,v.inRight)=rhs(end,v.inRight)-wall(2,v.inRight).*v.vx(v.inRight)';
H=zeros(size(g.H));
for j=1:numel(v.normal)
    speed=v.normal(j); ids=v.group{j};
    A=kron(speye(N),base.delta*M-speed*G+diag([max(-speed,0),max(speed,0)]));
    if speed>0
        A=A+sparse(3:2:2*N,2:2:2*N-2,-speed*ones(1,N-1),2*N,2*N);
    elseif speed<0
        A=A+sparse(2:2:2*N-2,3:2:2*N,speed*ones(1,N-1),2*N,2*N);
    end
    H(:,ids)=A\rhs(:,ids);
end
H=H-mean(H*v.weight);
assert(max(abs(H-r.H),[],'all')<1e-12);
fprintf('PASS second kinetic step uses the mixed macro source, not overwritten kinetic moments.\n');

% Evaluate the actual coupled DG residual independently from local weak
% transport, the complete numerical trace and the NEW source moments.
gain=(r.synthetic.sourceMoments/v.gram)*v.basis';
R=zeros(size(r.H));
for i=1:N
    ids=2*i-1:2*i;
    R(ids,:)=(-G*r.H(ids,:)).*v.vx' ...
        +[-r.kineticTrace(i,:);r.kineticTrace(i+1,:)].*v.vx' ...
        +M*(base.delta*(r.H(ids,:)-gain(ids,:))-2*r.config.acceleration*v.vy');
end
normR=sqrt(sum(sum(R.*(kron(speye(N),M)\R),1)'.*v.weight));
assert(abs(normR-r.diagnostics.coupled)<1e-12);
fprintf('PASS coupled residual vs independent kinetic weak residual.\n');

% Uniform translating equilibrium must remain an exact coupled fixed point.
eq=base; eq.flow='Couette'; eq.maxIter=10; eq.normalQuadrature='hermite';
eq.velocityOrder=[8 8]; eq.wallVelocity=[.007 .007];
eq.initial=@(x,vx,vy) 0*x+0*vx+2*.007*vy;
equilibrium=monoGSIS1D2V(eq);
assert(equilibrium.converged && equilibrium.iterations==1);
assert(max(abs(equilibrium.H-2*.007*equilibrium.velocity.vy'),[],'all')<1e-12);
fprintf('PASS uniform translating equilibrium, including the coupled source.\n');

% A tiny moment update alone must not imply a solved coupled equation.
short=base; short.delta=100; short.tol=1e-4;
notDone=monoGSIS1D2V(short);
assert(~notDone.converged && notDone.diagnostics.coupled>short.tol);
fprintf('PASS coupled residual prevents false stopping.\n');
fprintf('All GSIS source-coupling checks passed.\n');
end
