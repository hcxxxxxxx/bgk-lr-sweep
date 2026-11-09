function check_channel_reference
%CHECK_CHANNEL_REFERENCE Same-model reference and quadrature convergence.
root = fileparts(fileparts(mfilename('fullpath'))); addpath(root);
diary(fullfile(root,'debug','channel_reference_validation.log'));
cleanup = onCleanup(@() diary('off'));
results = mono1D2V(struct('plot',false,'save',false));
for name = {'couette','poiseuille'}
    r = results.(name{1}); cfg = r.config;
    ref = r.reference;
    coarse = channel1D2VReference(cfg,512,160);
    fineVelocity = channel1D2VReference(cfg,1024,320);
    finerSpace = channel1D2VReference(cfg,2048,160);
    spaceChange = norm(mean(reshape(ref.velocity,2,[]),1)'-coarse.velocity)/norm(coarse.velocity);
    velocityChange = norm(ref.velocity-fineVelocity.velocity)/norm(fineVelocity.velocity);
    finerSpaceChange = norm(mean(reshape(finerSpace.velocity,2,[]),1)'-ref.velocity)/norm(ref.velocity);
    fprintf('%s reference refinement: space 512->1024 %.5e, 1024->2048 %.5e; velocity 160->320 %.5e\n', ...
        name{1},spaceChange,finerSpaceChange,velocityChange);
    assert(spaceChange<2e-6 && finerSpaceChange<spaceChange && velocityChange<2e-6);
    assert(ref.linearResidual<1e-12);
    assert(r.converged && r.reference.relativeL2<1e-4);
    assert(max(r.diagnostics.balance)<1e-8);
    assert(max(abs(r.velocity.gram-diag([1 .5 .5 1])),[],'all')<1e-11);
    assert(max(abs(r.cellMoments(:,[1 2 4])),[],'all')<1e-12);
    % Deliberately reproduce the original full-range Hermite discretization.
    oldCfg = cfg; oldCfg.normalQuadrature = 'hermite';
    oldCfg.velocityOrder = [32 8];
    old = mono1D2V(oldCfg);
    assert(old.reference.relativeL2>0.02);
    % Refine space and velocity together: their errors can cancel, so
    % velocity-only refinement need not monotonically reduce the total error.
    refinedCfg = cfg; refinedCfg.velocityOrder = [256 8]; refinedCfg.cells = 128;
    refined = mono1D2V(refinedCfg);
    assert(refined.converged && refined.reference.relativeL2<r.reference.relativeL2);
    fprintf('PASS %s: old Hermite32/N64 %.6e, halfLegendre128/N64 %.6e, halfLegendre256/N128 %.6e (relative L2).\n', ...
        name{1},old.reference.relativeL2,r.reference.relativeL2,refined.reference.relativeL2);
    fprintf('Flow rate: default %.12g, reference %.12g; NS no slip %.12g.\n', ...
        r.flowRate,ref.flowRate,cfg.gap*mean(cfg.wallVelocity)+cfg.acceleration*cfg.delta*cfg.gap^3/6);
end
fprintf('All same-model BGK reference checks passed.\n');
end
