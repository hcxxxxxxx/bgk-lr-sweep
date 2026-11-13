function summary = check_reference_closure_1d2v
%CHECK_REFERENCE_CLOSURE_1D2V Isolate HOT and wall choices at a CIS fixed point.
% All errors below compare ONE macro solve with the same frozen CIS state.
% These are not errors from a converged, modified GSIS iteration.
% No production solver or stabilization parameter is changed.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
data=load(fullfile(root,'results','cubic_stopping_residuals32','comparison.mat'),'runs');
rows=cell(0,11);
for j=1:numel(data.runs)
    if data.runs(j).tolerance~=1e-9; continue; end
    for flow={'couette','poiseuille'}
        k=data.runs(j).cis.(flow{1}); v=k.velocity; N=k.mesh.count;
        delta=k.config.delta;
        [original,~,hdg,hot]=gsis_macro_from_kinetic_1d2v(k);
        assert(k.converged && isequal(hdg.cache.stabilization,[1 delta delta delta]));
        faceHOT=compute_high_order_moments_1d2v(k.kineticTrace,v);
        weak=project_hot_derivatives_1d(k.mesh,hot.R,faceHOT.R,'weak');
        interiorBC=macro_bc_from_kinetic(k.H([1 end],:),v);
        faceQ=k.kineticTrace*v.moment;
        completeBC=struct('trace',[faceQ([1 end],1)+faceQ([1 end],4), ...
            zeros(2,1),faceQ([1 end],3:4)]);
        param=struct('cache',hdg.cache);
        weakOnly=hdg_solve_macro_1d2v(k.mesh,delta,[0 k.config.acceleration], ...
            weak(:,1),weak(:,2),weak(:,3),interiorBC,param);
        wallOnly=hdg_solve_macro_1d2v(k.mesh,delta,[0 k.config.acceleration], ...
            hot.Pi_xx,hot.Pi_xy,hot.Theta_x,completeBC,param);
        weakAndWall=hdg_solve_macro_1d2v(k.mesh,delta,[0 k.config.acceleration], ...
            weak(:,1),weak(:,2),weak(:,3),completeBC,param);
        denominator=norm(k.cellMoments(:,3));
        errorPercent=@(m) 100*norm(m.cellAverage(:,3)-k.cellMoments(:,3))/denominator;

        % Insert kinetic moments and stresses into the macro weak equations.
        % The kinetic-derived trace is an illustrative embedding, not a claim
        % that the freely solved HDG trace must equal this particular trace.
        sigma=2*k.H*(v.weight.*v.vx.*v.vy);
        faceSigma=2*k.kineticTrace*(v.weight.*v.vx.*v.vy);
        uy=k.nodalMoments(:,3); faceUy=faceQ(:,3);
        stabilization=hdg.cache.stabilization(3);
        fluxDefect=zeros(2,N); auxiliaryWeak=0; auxiliaryStrong=0;
        for cellId=1:N
            ids=2*cellId-1:2*cellId;
            M=diff(k.mesh.edges(cellId:cellId+1))*[2 1;1 2]/6;
            G=[-1 -1;1 1]/2;
            fluxDefect(:,cellId)=[ ...
                -sigma(ids(1))+stabilization*(uy(ids(1))-faceUy(cellId))+faceSigma(cellId); ...
                 sigma(ids(2))+stabilization*(uy(ids(2))-faceUy(cellId+1))-faceSigma(cellId+1)];
            common=-delta*M*sigma(ids)+G*uy(ids)+[faceUy(cellId);-faceUy(cellId+1)];
            auxiliaryWeak=max(auxiliaryWeak,max(abs(common-M*weak(ids,2))));
            auxiliaryStrong=max(auxiliaryStrong,max(abs(common-M*hot.Pi_xy(ids))));
        end
        internal=max(abs([fluxDefect(2,1:end-1),fluxDefect(1,2:end)]));
        wall=max(abs([fluxDefect(1,1),fluxDefect(2,end)]));
        rows(end+1,:)={flow{1},delta,errorPercent(original),errorPercent(weakOnly), ...
            errorPercent(wallOnly),errorPercent(weakAndWall),internal,wall, ...
            auxiliaryStrong,auxiliaryWeak,mat2str(hdg.cache.stabilization)}; %#ok<AGROW>
    end
end
summary=cell2table(rows,'VariableNames',{'flow','delta','originalOneMacroUyPercent', ...
    'weakHOTOnlyPercent','completeWallOnlyPercent','weakHOTAndCompleteWallPercent', ...
    'internalFluxDefectWithKineticMomentTrace','wallFluxDefectWithKineticMomentTrace', ...
    'strongHOTAuxiliaryDefectWithKineticTrace','weakHOTAuxiliaryDefectWithKineticTrace', ...
    'stabilization'});
assert(height(summary)==6 && max(summary.weakHOTAuxiliaryDefectWithKineticTrace)<1e-12);
folder=fullfile(root,'results','reference_closure_comparison');
if ~exist(folder,'dir'); mkdir(folder); end
writetable(summary,fullfile(folder,'summary.csv'));
save(fullfile(folder,'diagnostics.mat'),'summary');
disp(summary);
fprintf('Auxiliary moment identity checked; full coupled fixed-point consistency remains unverified.\n');
end
