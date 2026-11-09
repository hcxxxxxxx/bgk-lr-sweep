function checks = check_delta01_reference_accuracy
% Isolate velocity, spatial and iteration errors against saved BGK reference.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
folder=fullfile(root,'results','delta01_reference_accuracy');
if ~exist(folder,'dir');mkdir(folder);end
diary(fullfile(folder,'run.log'));clean=onCleanup(@()diary('off')); %#ok<NASGU>
saved=load(fullfile(root,'results','channel_presentation_20261109','presentation_data.mat'),'cases');
rows=cell(0,11);details=struct([]);
% First isolate stopping, then normal velocity, then spatial resolution.
settings=[64 32 1e-5;64 32 1e-9;64 64 1e-9;64 128 1e-9;64 256 1e-9;128 256 1e-9;256 256 1e-9];
for flow={'couette','poiseuille'}
    id=find([saved.cases.delta]==.1 & [saved.cases.cells]==64 & strcmp({saved.cases.flow},flow{1}));
    base=saved.cases(id);cfg=base.cis.config;
    ref=channel1D2VReference(cfg,2048,320);
    for k=1:size(settings,1)
        N=settings(k,1);nv=settings(k,2);tol=settings(k,3);
        cfg.cells=N;cfg.velocityOrder=[nv 32];cfg.tol=tol;
        cfg.reference=false;cfg.plot=false;cfg.save=false;cfg.reportEvery=1000;
        if k==1;r=base.cis;else;r=mono1D2V(cfg);end
        assert(r.converged);
        uref=mean(reshape(ref.velocity,[],N),1)';u=r.cellMoments(:,3);
        scale=norm(uref);err=100*norm(u-uref)/scale;
        utarget=mean(reshape(r.cellMoments(:,3),[],64),1)';
        oldtarget=base.cis.cellMoments(:,3);
        qp=r.flowRate/(cfg.acceleration*cfg.gap^2);
        if cfg.acceleration==0;qp=NaN;end
        row={flow{1},N,nv,32,tol,r.iterations,err,100*norm(utarget-oldtarget)/norm(oldtarget),qp, ...
            min(abs(r.velocity.normal)),r.diagnostics.kinetic};
        rows(end+1,:)=row; %#ok<AGROW>
        details(end+1).result=r; %#ok<AGROW>
        fprintf('CHECK %s Nx=%d Nv=%dx32 tol=%.1e: reference error %.9g%%, Qp=%.12g\n', ...
            flow{1},N,nv,tol,err,qp);
    end
end
checks=cell2table(rows,'VariableNames',{'flow','cells','Nvx','Nvy','tolerance','iterations', ...
    'relativeVelocityReferencePercent','relativeChangeFromOriginalPercent','Qp', ...
    'minimumAbsVx','originalKineticResidual'});
writetable(checks,fullfile(folder,'summary.csv'));
save(fullfile(folder,'checks.mat'),'checks','details');
disp(checks);
end
