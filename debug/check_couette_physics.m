function check_couette_physics
%CHECK_COUETTE_PHYSICS Verify the requested Couette physical observables.
root = fileparts(fileparts(mfilename('fullpath'))); addpath(root);
diary(fullfile(root,'debug','couette_physics_validation.log'));
cleanup = onCleanup(@() diary('off'));
data = load(fullfile(root,'results','channel_1d2v','couette','channel.mat'),'r');
for delta = [1 10]
    if delta==1
        r = data.r;
    else
        % A separate, more collisional case; keep the delta=1 result intact.
        r = mono1D2V(struct('flow','Couette','delta',delta,'cells',128, ...
            'plot',false,'save',true,'outputDir',fullfile(root,'results','channel_1d2v','delta10')));
    end
    assert(r.converged && r.reference.relativeL2<2e-4);
    d = plotCouette1D2V(r);
    assert(d.maxNormalVelocity<1e-12);
    assert(d.maxDensityPerturbation<1e-12 && d.maxTemperaturePerturbation<1e-12);
    assert(d.shearStressRelativeSpan<1e-5);
    assert(d.actualWallSlip(1)>0 && d.actualWallSlip(2)<0);
    assert(max(abs(d.actualWallSlip+fliplr(d.actualWallSlip)))<1e-12);
    % Resolve curvature beyond the estimated DG/reference difference.
    assert(max(abs(d.fitRemainder))>10*r.reference.maximumAbsoluteError);
    assert(d.centralFitRelativeL2<0.01);
    assert(all(diff(r.cellMoments(:,3))>0));
    fprintf('PASS delta=%g: physical checks, slip and resolved curvature.\n',delta);
    close(gcf);
end
fprintf('All Couette physical checks passed.\n');
end
