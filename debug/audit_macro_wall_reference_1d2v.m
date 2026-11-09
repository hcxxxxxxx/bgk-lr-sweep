function summary = audit_macro_wall_reference_1d2v
%AUDIT_MACRO_WALL_REFERENCE_1D2V Check reference wall data and boundary sensitivity.
% Uses saved converged-CIS/one-HDG cases; does not change production BCs.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
folder=fullfile(root,'results','gsis_macro_1d2v','converged_cis_one_step');
data=load(fullfile(folder,'one_step.mat'),'runs');
diary(fullfile(folder,'wall_audit.log')); cleanup=onCleanup(@() diary('off'));
rows=cell(0,12);
for i=1:numel(data.runs)
    for name={'couette','poiseuille'}
        item=data.runs(i).(name{1}); k=item.kinetic;
        v=k.velocity; Hwall=k.H([1 end],:); w=v.weight;
        % Independent 2V adaptation of Fortran lines 538-550:
        rho=Hwall*w; temperature=Hwall*(w.*(v.vx.^2+v.vy.^2-1));
        uy=Hwall*(w.*v.vy);
        reference=[rho+temperature,zeros(2,1),uy,temperature];
        dataError=max(abs(reference-item.hot.boundary.trace),[],'all');
        appliedError=max(abs(reference-item.hdg.Uhat([1 end],:)),[],'all');
        difference=item.macro.uy-k.nodalMoments(:,3);
        relative=@(z) norm((z(1:2:end)+z(2:2:end))/2)/norm(k.cellMoments(:,3));
        % A change in wall uy data alone produces an affine homogeneous
        % response. Remove exactly the endpoint mismatch and measure what
        % remains. This is a diagnostic, NOT a proposed boundary condition.
        x=k.mesh.nodes; edges=k.mesh.edges;
        t=(x-edges(1))/(edges(end)-edges(1));
        endpointPart=(1-t)*difference(1)+t*difference(end);
        matched=difference-endpointPart;
        % Hold all volume terms fixed and only replace interior-wall
        % moments by complete diffuse/upwind numerical boundary moments.
        full=k.kineticTrace([1 end],:)*v.moment;
        fullBC=struct('trace',[full(:,1)+full(:,4),zeros(2,1),full(:,3:4)]);
        hot=item.hot;
        alternative=hdg_solve_macro_1d2v(k.mesh,k.config.delta, ...
            [0 k.config.acceleration],hot.Pi_xx,hot.Pi_xy,hot.Theta_x,fullBC, ...
            struct('cache',item.hdg.cache));
        fullChange=relative(alternative.uy-k.nodalMoments(:,3));
        assert(dataError<1e-12 && appliedError<1e-12);
        rows(end+1,:)={name{1},k.config.delta,dataError,appliedError, ...
            reference(1,3),reference(2,3),difference(1),difference(end), ...
            relative(difference),relative(matched),fullChange, ...
            max(abs(fullBC.trace-reference),[],'all')}; %#ok<AGROW>
        fprintf('%s delta=%g: reference wall uy=[%.12g %.12g], wall speeds=[%.12g %.12g]\n', ...
            name{1},k.config.delta,reference(:,3),k.config.wallVelocity);
        fprintf('  data mismatch %.3e, imposed trace mismatch %.3e; local endpoint uy difference=[%.3e %.3e]\n', ...
            dataError,appliedError,difference([1 end]));
        fprintf('  raw cell uy relative difference %.6e, after exact endpoint matching %.6e, complete-trace BC alternative %.6e\n', ...
            relative(difference),relative(matched),fullChange);
    end
end
summary=cell2table(rows,'VariableNames',{'flow','delta','referenceBoundaryDataError', ...
    'appliedBoundaryTraceError','leftGasUy','rightGasUy','leftLocalUyDifference', ...
    'rightLocalUyDifference','rawRelativeCellUyDifference', ...
    'relativeDifferenceAfterEndpointMatching','relativeDifferenceWithCompleteTraceBC', ...
    'maximumCompleteVsInteriorBoundaryDifference'});
writetable(summary,fullfile(folder,'wall_audit.csv'));
disp(summary);
end
