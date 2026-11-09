function report = audit_hdg_weak_form_1d2v
%AUDIT_HDG_WEAK_FORM_1D2V Audit operators against independently integrated PDEs.
% Does not construct the reference operator using production F/G/M, kron,
% response matrices or cached face-flux matrices. Tests each basis vector
% against scalar weak equations and explicit numerical face flux formulas.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
folder=fullfile(root,'results','gsis_macro_1d2v');
if ~exist(folder,'dir'); mkdir(folder); end
diary(fullfile(root,'debug','hdg_weak_form_audit.log'));
cleanup=onCleanup(@() diary('off'));
mesh=struct('edges',[-.5;-.29;.03;.5]);
strengths=[.07 1 7 1000]; errors=zeros(8,8); row=0;
for delta=strengths
    settings={[1 delta delta delta],[.7 1.1 .4 .9]};
    for j=1:numel(settings)
        stabilization=settings{j}; param=struct('stabilization',stabilization);
        if j==1; param=struct; end % Audit the actual default as well as an override.
        cache=hdg_precompute_p1(mesh,delta,param);
        assert(isequal(cache.stabilization,stabilization));
        matrixError=0; fluxError=0; schurError=0; sourceError=0; strongBlockError=0;
        for cellId=1:cache.count
            lengthCell=diff(mesh.edges(cellId:cellId+1));
            [A,B,C,D]=reference_operators(lengthCell,delta,stabilization);
            matrixError=max(matrixError,max(abs([A-cache.A(:,:,cellId), ...
                B-cache.B(:,:,cellId)]),[],'all'));
            fluxError=max(fluxError,max(abs([C-cache.C(:,:,cellId), ...
                D-cache.D(:,:,cellId)]),[],'all'));
            Kref=C*(A\B)+D;
            Kcode=cache.C(:,:,cellId)*cache.response(:,:,cellId)+cache.D(:,:,cellId);
            schurError=max(schurError,max(abs(Kref-Kcode),[],'all'));
            forcing=reshape(sin(1:14),2,7);
            loadReference=-reference_weak(zeros(14,1),zeros(8,1),lengthCell, ...
                delta,stabilization,forcing);
            loadCode=reshape(cache.M(:,:,cellId)*forcing,[],1);
            sourceError=max(sourceError,max(abs(loadReference-loadCode)));
            % Integrating the auxiliary face term back into the volume must
            % give the strong derivative block, as in S-model Eq.(24).
            Cphi=[-1 1;-1 1]/2; expected=zeros(8,6);
            expected(3:4,1:2)=-Cphi/delta;
            expected(5:6,3:4)=-Cphi/delta;
            expected(7:8,5:6)=-Cphi/delta;
            strongBlockError=max(strongBlockError, ...
                max(abs(cache.A(1:8,9:14,cellId)-expected),[],'all'));
        end
        row=row+1; errors(row,:)=[delta,stabilization(1),matrixError,fluxError, ...
            schurError,sourceError,min(cache.localRcond),strongBlockError];
        assert(matrixError<1e-12 && fluxError<1e-12 && sourceError<1e-12);
        assert(schurError<1e-10 && strongBlockError<1e-12);
    end
end
operatorAudit=array2table(errors,'VariableNames',{'delta','sRho','localABError', ...
    'faceCDError','localSchurError','sourceError','minimumLocalRcond','auxiliaryStrongBlockError'});
writetable(operatorAudit,fullfile(folder,'operator_audit.csv'));
fprintf('Independent integration: max AB %.3e, CD %.3e, Schur %.3e, source %.3e\n', ...
    max(errors(:,3:6),[],1));

% Check the discrete global recovery and boundary sign using new reference
% operators and independently evaluated weak/flux residuals.
param=struct('stabilization',[.7 1.1 .4 .9]); delta=2.3;
x=reshape([mesh.edges(1:end-1),mesh.edges(2:end)]',[],1);
acceleration=[.013+.004*x,-.007+.003*x];
correction=[.02+.03*x,.01-.04*x,.003+.005*x];
bc=struct('trace',[.01 .002 -.01 .003;.02 .002 NaN .006], ...
    'flux',[NaN NaN NaN NaN;NaN NaN .004 NaN]);
[~,~,hdg]=hdg_solve_macro_1d2v(mesh,delta,acceleration, ...
    correction(:,1),correction(:,2),correction(:,3),bc,param);
forcing=[zeros(numel(x),1),2*acceleration,zeros(numel(x),1),correction];
globalReference=zeros(4*(hdg.cache.count+1)); rhsReference=zeros(size(globalReference,1),1);
independentLocalResidual=0;
for i=1:hdg.cache.count
    h=diff(mesh.edges(i:i+1)); ids=hdg.cache.traceIds(:,i);
    [A,B,C,D]=reference_operators(h,delta,param.stabilization);
    load=-reference_weak(zeros(14,1),zeros(8,1),h,delta,param.stabilization,forcing(2*i-1:2*i,:));
    globalReference(ids,ids)=globalReference(ids,ids)+C*(A\B)+D;
    rhsReference(ids)=rhsReference(ids)-C*(A\load);
    local=reference_weak(hdg.localFields(:,i),reshape(hdg.Uhat(i:i+1,:)',[],1), ...
        h,delta,param.stabilization,forcing(2*i-1:2*i,:));
    independentLocalResidual=max(independentLocalResidual,max(abs(local)));
end
globalMatrixError=max(abs(globalReference-full(hdg.cache.K)),[],'all');
globalRhsError=max(abs(rhsReference-hdg.traceRhs));
assert(globalMatrixError<1e-12 && globalRhsError<1e-12 && independentLocalResidual<1e-12);
% Recover the explicit interface-mean formulas from flux continuity. This
% is the 1D2V counterpart of Eq.(30) in the original S-model note.
traceFormulaError=0; U=hdg.localFields; s=param.stabilization;
for i=1:hdg.cache.count-1
    expected=(U([2 4 6 8],i)+U([1 3 5 7],i+1))/2;
    expected(1)=expected(1)+(U(4,i)-U(3,i+1))/(2*s(1));
    expected(2)=expected(2)-(U(10,i)-U(9,i+1))/(2*delta*s(2));
    expected(3)=expected(3)-(U(12,i)-U(11,i+1))/(2*delta*s(3));
    expected(4)=expected(4)-(U(14,i)-U(13,i+1))/(2*delta*s(4));
    traceFormulaError=max(traceFormulaError,max(abs(expected-hdg.Uhat(i+1,:)')));
end
boundaryResidual=0;
for side=1:2
    if side==1; i=1; else; i=hdg.cache.count; end
    flux=reshape(reference_face_flux(hdg.localFields(:,i), ...
        reshape(hdg.Uhat(i:i+1,:)',[],1),delta,s),4,2);
    node=1+(side-1)*hdg.cache.count;
    for component=1:4
        if isfinite(bc.trace(side,component))
            value=hdg.Uhat(node,component)-bc.trace(side,component);
        else
            value=flux(component,side)-bc.flux(side,component);
        end
        boundaryResidual=max(boundaryResidual,abs(value));
    end
end
assert(traceFormulaError<1e-12 && boundaryResidual<1e-12);
fprintf('Independent global assembly: K %.3e, rhs %.3e, weak residual %.3e\n', ...
    globalMatrixError,globalRhsError,independentLocalResidual);
fprintf('Interface formulas %.3e; boundary Dirichlet/flux residual %.3e; strong auxiliary blocks %.3e\n', ...
    traceFormulaError,boundaryResidual,max(errors(:,8)));

% Verify 2V Gaussian projection/stress/heat coefficients and force moments.
velocity=reference_velocity;
v=velocity; w=v.weight; r2=v.vx.^2+v.vy.^2;
invariant=[ones(numel(w),1),v.vx,v.vy,r2-1];
stressHeat=[v.vx.^2-v.vy.^2,2*v.vx.*v.vy,v.vx.*(r2-2)];
Gram=invariant'*(w.*invariant);
projectionError=max(abs(Gram-diag([1 .5 .5 1])),[],'all');
orthogonality=max(abs(stressHeat'*(w.*invariant)),[],'all');
H=[1,2,-3,.7]*invariant'+.13*(v.vx.*v.vy)';
Q=H*(w.*invariant); moments=H*(w.*stressHeat);
conservationIdentity=H*(w.*v.vx.*invariant)- ...
    [Q(2),(Q(1)+Q(4)+moments(1))/2,moments(2)/2,Q(2)+moments(3)];
hot=compute_high_order_moments_1d2v(H,v);
fluxIdentity=[H*(w.*v.vx.*stressHeat(:,1)), ...
    H*(w.*v.vx.*stressHeat(:,2)),H*(w.*v.vx.*stressHeat(:,3))] ...
    -[Q(2)+hot.Rxx,Q(3)+hot.Rxy,Q(4)+hot.RT];
force=2*(.031*v.vx-.027*v.vy);
forceMoments=force'*(w.*invariant);
forceError=max(abs(forceMoments-[0 .031 -.027 0]));
assert(max([projectionError,orthogonality,max(abs(fluxIdentity)), ...
    max(abs(conservationIdentity)),forceError])<1e-12);
fprintf('Gaussian coefficients: Gram %.3e, stress/heat orthogonality %.3e, HOT identity %.3e, force %.3e\n', ...
    projectionError,orthogonality,max(abs(fluxIdentity)),forceError);

% Locate a kinetic/HDG fixed-point mismatch BEFORE solving the macro problem.
% Insert kinetic moments/stresses, illustrative kinetic-derived interior
% traces and reference-code wall data. This is a diagnostic, not a demand
% that the free HDG interior traces equal kinetic upwind moments.
cases=mono1D2V(struct('plot',false,'save',false,'reference',false, ...
    'stopCriterion','equation','tol',1e-10));
fixedPoint=zeros(2,10); names={'couette','poiseuille'};
for j=1:2
    k=cases.(names{j}); v=k.velocity; w=v.weight; N=k.mesh.count;
    [~,~,hdg,hot]=gsis_macro_from_kinetic_1d2v(k);
    B=[ones(numel(w),1),2*v.vx,2*v.vy,v.speed2-2];
    localStress=k.H*(w.*[v.vx.^2-v.vy.^2,2*v.vx.*v.vy,v.vx.*(v.speed2-2)]);
    fluxLocal=k.H*(w.*v.vx.*B);
    hat=k.kineticTrace*v.moment; hat(:,1)=hat(:,1)+hat(:,4);
    hat([1 end],:)=hot.boundary.trace;
    fluxKinetic=k.kineticTrace*(w.*v.vx.*B);
    forceNodal=2*k.config.acceleration*v.vy';
    forceMoment=forceNodal*(w.*B);
    primitive=k.nodalMoments; primitive(:,1)=primitive(:,1)+primitive(:,4);
    fields=[primitive,-k.config.delta*localStress];
    s=hdg.cache.stabilization; kineticResidual=0; conservationResidual=0; auxiliaryResidual=0;
    fluxDifference=0; decompositionError=0; boundaryMatch=0;
    for i=1:N
        ids=2*i-1:2*i; h=diff(k.mesh.edges(i:i+1));
        U=reshape(fields(ids,:),[],1); hatCell=reshape(hat(i:i+1,:)',[],1);
        forcing=zeros(2,7); forcing(:,3)=2*k.config.acceleration;
        forcing(:,5:7)=hot.derivative(ids,:);
        residual=reshape(reference_weak(U,hatCell,h,k.config.delta,s,forcing),2,7);
        conservationResidual=max(conservationResidual,max(abs(residual(:,1:4)),[],'all'));
        auxiliaryResidual=max(auxiliaryResidual,max(abs(residual(:,5:7)),[],'all'));
        M=h*[2 1;1 2]/6; G=[-1 -1;1 1]/2;
        kinetic=-G*fluxLocal(ids,:)+[-fluxKinetic(i,:);fluxKinetic(i+1,:)] ...
            -M*repmat(forceMoment,2,1);
        kineticResidual=max(kineticResidual,max(abs(kinetic),[],'all'));
        fluxHDG=reference_face_flux(U,hatCell,k.config.delta,s);
        fluxOther=[-fluxKinetic(i,:)';fluxKinetic(i+1,:)'];
        difference=reshape(fluxHDG-fluxOther,4,2)';
        fluxDifference=max(fluxDifference,max(abs(difference),[],'all'));
        sourceDifference=M*(repmat(forceMoment,2,1)-forcing(:,1:4));
        check=residual(:,1:4)-kinetic-difference-sourceDifference;
        decompositionError=max(decompositionError,max(abs(check),[],'all'));
        if i==1; boundaryMatch=max(boundaryMatch,max(abs(hatCell(1:4)-hot.boundary.trace(1,:)'))); end
        if i==N; boundaryMatch=max(boundaryMatch,max(abs(hatCell(5:8)-hot.boundary.trace(2,:)'))); end
    end
    fixedPoint(j,:)=[k.config.delta,kineticResidual,conservationResidual,auxiliaryResidual, ...
        fluxDifference,decompositionError,boundaryMatch,max(abs(hot.velocityMomentDefect)), ...
        k.diagnostics.kinetic,max(abs(forceMoment-[0 0 2*k.config.acceleration 0]))];
    assert(decompositionError<1e-12 && boundaryMatch<1e-12);
    fprintf('%s fixed-point probe: kinetic moment %.3e; HDG conservation %.3e; auxiliary %.3e; face difference %.3e; decomposition %.3e\n', ...
        names{j},kineticResidual,conservationResidual,auxiliaryResidual,fluxDifference,decompositionError);
end
fixedPointAudit=array2table(fixedPoint,'VariableNames',{'delta','kineticMomentWeakResidual', ...
    'hdgConservationWeakResidual','hdgAuxiliaryWeakResidual','maximumFaceFluxDifference', ...
    'residualDecompositionError','boundaryTraceMismatch','velocityMomentDefect', ...
    'kineticResidual','discreteForceMomentDefect'});
fixedPointAudit.flow=names'; fixedPointAudit=movevars(fixedPointAudit,'flow','Before',1);
writetable(fixedPointAudit,fullfile(folder,'fixed_point_audit.csv'));
report=struct('operators',operatorAudit,'fixedPoint',fixedPointAudit, ...
    'globalMatrixError',globalMatrixError,'globalRhsError',globalRhsError, ...
    'traceFormulaError',traceFormulaError,'boundaryResidual',boundaryResidual, ...
    'independentWeakResidual',independentLocalResidual,'projectionError',projectionError, ...
    'orthogonalityError',orthogonality,'forceMomentError',forceError, ...
    'conservationMomentIdentityError',max(abs(conservationIdentity)), ...
    'highOrderMomentIdentityError',max(abs(fluxIdentity)));
save(fullfile(folder,'weak_form_audit.mat'),'report');
fprintf('HDG formula audit passed. Fixed-point differences are diagnostics, not a passed equivalence test.\n');
end

function [A,B,C,D]=reference_operators(h,delta,s)
A=zeros(14); B=zeros(14,8); C=zeros(8,14); D=zeros(8);
for j=1:14
    U=zeros(14,1); U(j)=1;
    A(:,j)=reference_weak(U,zeros(8,1),h,delta,s,zeros(2,7));
    C(:,j)=reference_face_flux(U,zeros(8,1),delta,s);
end
for j=1:8
    hat=zeros(8,1); hat(j)=1;
    B(:,j)=-reference_weak(zeros(14,1),hat,h,delta,s,zeros(2,7));
    D(:,j)=reference_face_flux(zeros(14,1),hat,delta,s);
end
end

function residual=reference_weak(U,hat,h,delta,s,forcing)
% Direct quadrature of the seven scalar weak equations.
field=reshape(U,2,7); boundary=reshape(hat,4,2);
points=[-sqrt(3/5),0,sqrt(3/5)]; weights=[5/9,8/9,5/9]*h/2;
derivative=[-1;1]/h; residual=zeros(2,7);
for p=1:3
    shape=[1-points(p);1+points(p)]/2; value=shape'*field;
    flux=[value(2),value(1)-value(5)/delta,-value(6)/delta,-value(7)/delta];
    residual(:,1:4)=residual(:,1:4)-weights(p)*derivative*flux;
    residual(:,5:7)=residual(:,5:7)+weights(p)*(shape*value(5:7)+derivative*value(2:4));
    residual=residual-weights(p)*shape*(shape'*forcing);
end
face=reshape(reference_face_flux(U,hat,delta,s),4,2);
for side=1:2
    normal=2*side-3;
    residual(side,1:4)=residual(side,1:4)+face(:,side)';
    residual(side,5:7)=residual(side,5:7)-normal*boundary(2:4,side)';
end
residual=residual(:);
end

function flux=reference_face_flux(U,hat,delta,s)
local=reshape(U,2,7); trace=reshape(hat,4,2); flux=zeros(4,2);
for side=1:2
    n=2*side-3; q=local(side,:); z=trace(:,side);
    flux(1,side)=q(2)*n+s(1)*(q(1)-z(1));
    flux(2,side)=(z(1)-q(5)/delta)*n+s(2)*(q(2)-z(2));
    flux(3,side)=-q(6)/delta*n+s(3)*(q(3)-z(3));
    flux(4,side)=-q(7)/delta*n+s(4)*(q(4)-z(4));
end
flux=flux(:);
end

function v=reference_velocity
% Gaussian quadrature exact for all polynomial products used in this audit.
offdiag=sqrt((1:11)/2); [V,D]=eig(diag(offdiag,1)+diag(offdiag,-1));
[node,ids]=sort(diag(D)); weight=V(1,ids).^2;
[vx,vy]=meshgrid(node,node); w=weight'*weight;
v=struct('vx',vx(:),'vy',vy(:),'weight',w(:));
end
