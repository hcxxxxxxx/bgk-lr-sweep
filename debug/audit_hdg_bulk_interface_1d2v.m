function summary = audit_hdg_bulk_interface_1d2v(folder)
%AUDIT_HDG_BULK_INTERFACE_1D2V Independent HOT/face checks and exact macro solve.
% Loads the converged-CIS -> one-HDG cases. No production solver changes.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
if nargin<1; folder=fullfile(root,'results','gsis_macro_1d2v','converged_cis_one_step'); end
data=load(fullfile(folder,'one_step.mat'),'runs');
diary(fullfile(folder,'bulk_interface_audit.log')); cleanup=onCleanup(@() diary('off'));
rows=cell(0,16); details=struct([]); entry=0;
for i=1:numel(data.runs)
    for name={'couette','poiseuille'}
        entry=entry+1; item=data.runs(i).(name{1}); k=item.kinetic;
        v=k.velocity; w=v.weight; H=k.H; mesh=k.mesh; N=mesh.count;
        delta=k.config.delta; a=k.config.acceleration; s=item.hdg.cache.stabilization;
        edges=mesh.edges; dx=diff(edges); x=mesh.nodes;
        uy=k.nodalMoments(:,3); sigma=2*H*(w.*v.vx.*v.vy);
        check_reference_split(delta,s(3));
        higher=H*(w.*(v.vy.*(2*v.vx.^2-1)));
        J=2*H*(w.*v.vx.^2.*v.vy);
        Jhat=2*k.kineticTrace*(w.*v.vx.^2.*v.vy);
        sigmaHat=2*k.kineticTrace*(w.*v.vx.*v.vy);
        gain=(k.nodalMoments/v.gram)*v.basis'; force=2*a*v.vy';
        piIndependent=zeros(2*N,1); hotLoadError=0; auxiliaryIdentityError=0;
        strongClosure=0; kineticShearWeak=0; closureDecompositionError=0;
        momentumStrong=0; momentumDecompositionError=0; faceResidual=0;
        for cellId=1:N
            ids=2*cellId-1:2*cellId; h=dx(cellId);
            M=h*[2 1;1 2]/6; G=[-1 -1;1 1]/2;
            % Integrate phi_m * dH/dx directly in velocity space, as Fortran
            % INT_NODFUNC_TRI_TRI_X(L,M) does, without the production R helper.
            gradient=(H(ids(2),:)-H(ids(1),:))/h;
            piValue=gradient*(w.*v.vy.*(2*v.vx.^2-1));
            piIndependent(ids)=piValue;
            hotLoad=[h/2;h/2]*piValue;
            hotLoadError=max(hotLoadError,max(abs(hotLoad-M*item.hot.Pi_xy(ids))));
            sourceShear=(delta*gain(ids,:)+force)*(2*w.*v.vx.*v.vy);
            rKin=-G*J(ids)+[-Jhat(cellId);Jhat(cellId+1)] ...
                +delta*M*sigma(ids)-M*sourceShear;
            rStrong=M*(delta*sigma(ids)+piValue)+G'*uy(ids);
            lift=[-(Jhat(cellId)-J(ids(1)));Jhat(cellId+1)-J(ids(2))];
            kineticShearWeak=max(kineticShearWeak,max(abs(rKin)));
            strongClosure=max(strongClosure,max(abs(rStrong)));
            closureDecompositionError=max(closureDecompositionError, ...
                max(abs(rStrong-rKin+lift-M*sourceShear)));
            % Independently check the scaled y-momentum strong/weak identity.
            rMomentum=-G*sigma(ids)+[-sigmaHat(cellId);sigmaHat(cellId+1)]-M*[2*a;2*a];
            rMomentumStrong=G'*sigma(ids)-M*[2*a;2*a];
            momentumLift=[-(sigmaHat(cellId)-sigma(ids(1)));sigmaHat(cellId+1)-sigma(ids(2))];
            momentumStrong=max(momentumStrong,max(abs(rMomentumStrong)));
            momentumDecompositionError=max(momentumDecompositionError, ...
                max(abs(rMomentumStrong-rMomentum+momentumLift)));
            % Verify actual recovered HDG auxiliary equations from the
            % scalar reference weak form, not the cached A/B blocks.
            U=item.hdg.localFields(:,cellId); L=U(11:12);
            z=item.hdg.Uhat(cellId:cellId+1,3);
            rAux=M*(L-piValue)+G*U(5:6)+[z(1);-z(2)];
            auxiliaryIdentityError=max(auxiliaryIdentityError,max(abs(rAux)));
        end
        % Independent reference interface formula from the sum of outward
        % stresses plus penalties: uhat=mean(u)-[L]/(2*delta*s_y).
        U=item.hdg.localFields;
        for face=1:N-1
            expected=(U(6,face)+U(5,face+1))/2 ...
                -(U(12,face)-U(11,face+1))/(2*delta*s(3));
            faceResidual=max(faceResidual,abs(expected-item.hdg.Uhat(face+1,3)));
        end
        % Solve the CONTINUOUS frozen-HOT macro problem exactly. Pi is
        % constant in each P1 kinetic cell: L'= -2*delta*a; u'=L-Pi.
        piCell=piIndependent(1:2:end);
        uLeft=item.hot.boundary.trace(1,3); uRight=item.hot.boundary.trace(2,3);
        gap=edges(end)-edges(1);
        Lleft=(uRight-uLeft+sum(dx.*piCell))/gap+delta*a*gap;
        cumulative=[0;cumsum(dx.*piCell)];
        distance=x-edges(1);
        integralPi=repelem(cumulative(1:end-1),2) ...
            +piIndependent.*(x-repelem(edges(1:end-1),2));
        exactNodal=uLeft+Lleft*distance-delta*a*distance.^2-integralPi;
        % Exact quadratic cell average, not the trapezoidal mean of endpoints.
        left=edges(1:end-1)-edges(1); center=left+dx/2;
        exactCell=uLeft+Lleft*center-delta*a*(left.^2+left.*dx+dx.^2/3) ...
            -cumulative(1:end-1)-piCell.*dx/2;
        denom=norm(k.cellMoments(:,3));
        rawDifference=norm(item.macro.cellAverage(:,3)-k.cellMoments(:,3))/denom;
        exactVsCIS=norm(exactCell-k.cellMoments(:,3))/denom;
        hdgVsExact=norm(item.macro.cellAverage(:,3)-exactCell)/denom;
        nodalVsExact=max(abs(item.macro.uy-exactNodal));
        assert(hotLoadError<1e-12 && auxiliaryIdentityError<1e-12 && faceResidual<1e-12);
        assert(closureDecompositionError<1e-12 && momentumDecompositionError<1e-12);
        assert(max(abs(J-uy-higher))<1e-12);
        rows(end+1,:)={name{1},delta,hotLoadError,auxiliaryIdentityError,faceResidual, ...
            kineticShearWeak,strongClosure,closureDecompositionError,momentumStrong, ...
            momentumDecompositionError,rawDifference,exactVsCIS,hdgVsExact,nodalVsExact, ...
            max(abs(diff(reshape(higher,2,N),1,1))),mat2str(s)}; %#ok<AGROW>
        details(entry).flow=name{1}; details(entry).delta=delta;
        details(entry).x=x; details(entry).exactNodalUy=exactNodal;
        details(entry).exactCellUy=exactCell; details(entry).piIndependent=piIndependent;
        fprintf('%s delta=%g: HOT load %.3e, auxiliary %.3e, interface %.3e\n', ...
            name{1},delta,hotLoadError,auxiliaryIdentityError,faceResidual);
        fprintf('  HDG/CIS %.6e, exact frozen-HOT macro/CIS %.6e, HDG/exact macro %.6e\n', ...
            rawDifference,exactVsCIS,hdgVsExact);
        fprintf('  kinetic shear weak residual %.3e, broken strong closure defect %.3e, face-lift identity error %.3e\n', ...
            kineticShearWeak,strongClosure,closureDecompositionError);
    end
end
summary=cell2table(rows,'VariableNames',{'flow','delta','hotLoadError','auxiliaryWeakError', ...
    'interfaceFormulaError','kineticShearWeakResidual','brokenStrongClosureDefect', ...
    'closureFaceLiftIdentityError','brokenStrongMomentumDefect','momentumFaceLiftIdentityError', ...
    'rawHDGVsCIS','exactFrozenHOTMacroVsCIS','HDGVsExactFrozenHOTMacro','maximumNodalHDGVsExact', ...
    'maximumCellRxyVariation','stabilization'});
writetable(summary,fullfile(folder,'bulk_interface_audit.csv'));
save(fullfile(folder,'bulk_interface_audit.mat'),'summary','details');
disp(summary);
end

function check_reference_split(delta,s)
% Fortran has TWO auxiliaries, each receiving Pi_xy/2 (lines 449-450).
% At d/dy=0, Lxy+Lyx=uy_x+Pi_xy: the reduced 1D source is FULL Pi_xy.
h=.17; M=h*[2 1;1 2]/6; G=[-1 -1;1 1]/2;
piValue=[.013;-.007]; uhat=[-.02;.03]; acceleration=.011;
A=[s*eye(2),-G'/delta,-G'/delta;zeros(2),M,zeros(2);G,zeros(2),M];
rhs=[M*[2*acceleration;2*acceleration]+s*uhat; ...
    M*piValue/2;M*piValue/2+[-uhat(1);uhat(2)]];
reference=A\rhs;
reduced=[s*eye(2),-G'/delta;G,M]\ ...
    [M*[2*acceleration;2*acceleration]+s*uhat;M*piValue+[-uhat(1);uhat(2)]];
assert(max(abs(reference(1:2)-reduced(1:2)))<1e-12);
assert(max(abs(reference(3:4)+reference(5:6)-reduced(3:4)))<1e-12);
end
