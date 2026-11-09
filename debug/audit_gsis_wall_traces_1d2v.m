function summary = audit_gsis_wall_traces_1d2v
%AUDIT_GSIS_WALL_TRACES_1D2V Audit saved, fully coupled GSIS wall traces.
% Independently reconstruct boundary moments, diffuse inflows, upwind traces
% and outward macro fluxes. Does not modify or rerun the flow solvers.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
folder=fullfile(root,'results','gsis_channel_comparison','reference_stopping');
data=load(fullfile(folder,'comparison.mat'),'runs');
old=load(fullfile(root,'results','gsis_macro_1d2v', ...
    'converged_cis_one_step','one_step.mat'),'runs');
diary(fullfile(folder,'wall_trace_audit.log'));
cleanup=onCleanup(@() diary('off'));
rows=cell(0,16);
for i=1:numel(data.runs)
    for name={'couette','poiseuille'}
        g=data.runs(i).gsis.(name{1}); c=data.runs(i).cis.(name{1});
        H=g.H; v=g.velocity; w=v.weight; vx=v.vx; vy=v.vy;
        N=g.mesh.count; delta=g.config.delta;
        hwall=H([1 end],:);
        % Reference Fortran projects the INTERIOR VDF onto the wall face,
        % integrates over ALL velocities, then removes normal velocity.
        rho=hwall*w; T=hwall*(w.*(vx.^2+vy.^2-1)); uy=hwall*(w.*vy);
        target=[rho+T,zeros(2,1),uy,T];
        hdg=g.synthetic.hdg; hat=hdg.Uhat([1 end],:);
        dataError=max(abs(target-g.synthetic.hot.boundary.trace),[],'all');
        appliedError=max(abs(target-hat),[],'all');
        known=[1:4,4*N+(1:4)]';
        assert(isequal(double(hdg.bc.known(:)),known));
        assert(isequal(hdg.bc.free(:),setdiff((1:4*(N+1))',known)));
        assert(isequal(g.config.wallVelocity,c.config.wallVelocity));
        assert(g.config.acceleration==c.config.acceleration);
        assert(isequal(hdg.cache.stabilization,[1 delta delta delta]));

        left=vx>0; right=vx<0;
        nw=[sum(H(1,right)'.*w(right).*abs(vx(right)))/sum(w(left).*vx(left)); ...
            sum(H(end,left)'.*w(left).*abs(vx(left)))/sum(w(right).*abs(vx(right)))];
        diffuse=nw+2*g.config.wallVelocity(:)*vy';
        wallTrace=g.kineticTrace([1 end],:);
        inflowError=max([abs(wallTrace(1,left)-diffuse(1,left)), ...
            abs(wallTrace(2,right)-diffuse(2,right))]);
        outflowError=max([abs(wallTrace(1,right)-H(1,right)), ...
            abs(wallTrace(2,left)-H(end,left))]);
        internal=zeros(N-1,numel(vx));
        internal(:,left)=H(2:2:end-2,left);
        internal(:,right)=H(3:2:end,right);
        upwindError=max(abs(internal-g.kineticTrace(2:end-1,:)),[],'all');
        massFlux=max(abs(wallTrace*(w.*vx)));

        local=g.synthetic.macro.pressureNodal([1 end],:);
        L=[g.synthetic.aux.Lx([1 end]),g.synthetic.aux.Ly([1 end]), ...
            g.synthetic.aux.E([1 end])];
        normals=[-1;1]; s=hdg.cache.stabilization;
        physical=[local(:,2),hat(:,1)-L(:,1)/delta, ...
            -L(:,2)/delta,-L(:,3)/delta];
        flux=normals.*physical+(local-hat).*s;
        savedFlux=[hdg.outwardFlux(1:4,1)';hdg.outwardFlux(5:8,end)'];
        fluxError=max(abs(flux-savedFlux),[],'all');
        shear=wallTrace*(w.*vx.*vy);
        wallShearGap=normals.*flux(:,3)/2-shear;
        fullUy=wallTrace*(w.*vy);

        index=find([old.runs.delta]==delta,1); assert(~isempty(index));
        previous=old.runs(index).(name{1}); k=previous.kinetic;
        assert(isequal(k.config.wallVelocity,g.config.wallVelocity));
        assert(k.config.acceleration==g.config.acceleration);
        oldWall=k.H([1 end],:); ov=k.velocity;
        oldRho=oldWall*ov.weight;
        oldT=oldWall*(ov.weight.*(ov.vx.^2+ov.vy.^2-1));
        oldTarget=[oldRho+oldT,zeros(2,1),oldWall*(ov.weight.*ov.vy),oldT];
        oldError=max(abs(oldTarget-previous.hdg.Uhat([1 end],:)),[],'all');
        assert(max([dataError,appliedError,inflowError,outflowError, ...
            upwindError,massFlux,fluxError,oldError])<1e-12);
        rows(end+1,:)={name{1},delta,dataError,appliedError,inflowError, ...
            outflowError,upwindError,massFlux,fluxError,oldError, ...
            uy(1),uy(2),max(abs(local(:,3)-hat(:,3))), ...
            max(abs(fullUy-uy)),wallShearGap(1),wallShearGap(2)}; %#ok<AGROW>
    end
end
summary=cell2table(rows,'VariableNames',{'flow','delta','macroBoundaryDataError', ...
    'macroAppliedTraceError','currentDiffuseInflowError','outflowCopyError', ...
    'interiorUpwindTraceError','wallMassFlux','macroWallFluxFormulaError', ...
    'previousBoundaryRuleError','leftGasUy','rightGasUy', ...
    'localVsHatUyWallDifference','completeVsInteriorUyWallDifference', ...
    'leftMacroMinusKineticShear','rightMacroMinusKineticShear'});
writetable(summary,fullfile(folder,'wall_trace_audit.csv'));
disp(summary);
fprintf('PASS current GSIS wall traces, reference boundary rule, and previous saved boundary rule.\n');
end
