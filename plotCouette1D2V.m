function d = plotCouette1D2V(r,saveOutput)
%PLOTCOUETTE1D2V Couette physics checks, central fit and wall zoom.
%   plotCouette1D2V;                % inspect the saved default Couette run
%   d = plotCouette1D2V(r,false);    % inspect a returned single-flow result
% The central fit uses |x|<=0.2*gap. Its remainder is a diagnostic, not an
% exact half-space Knudsen-layer function: at high rarefaction the two wall
% layers can overlap, leaving no well-separated hydrodynamic bulk region.
if nargin<1
    root = fileparts(mfilename('fullpath'));
    data = load(fullfile(root,'results','channel_1d2v','couette','channel.mat'),'r');
    r = data.r;
end
if nargin<2; saveOutput = true; end
assert(strcmpi(r.config.flow,'Couette') && r.config.acceleration==0, ...
    'This diagnostic requires pure Couette flow with zero acceleration.');
assert(r.mesh.count>=4,'Use at least four cells for the central linear fit.');
cfg = r.config; x = r.mesh.center; q = r.cellMoments; u = q(:,3);
wallX = [-cfg.gap/2;cfg.gap/2]; gasU = r.nodalMoments([1 end],3);
central = abs(x)<=0.2*cfg.gap;
coefficient = [ones(nnz(central),1),x(central)]\u(central);
fit = coefficient(1)+coefficient(2)*x;
wallFit = coefficient(1)+coefficient(2)*wallX;
d = struct('maxNormalVelocity',max(abs(r.nodalMoments(:,2))), ...
    'maxDensityPerturbation',max(abs(r.nodalMoments(:,1))), ...
    'maxTemperaturePerturbation',max(abs(r.nodalMoments(:,4))), ...
    'meanShearStress',mean(r.shearStress), ...
    'shearStressSpan',max(r.shearStress)-min(r.shearStress), ...
    'centralFitInterval',[-0.2 0.2]*cfg.gap, ...
    'fitIntercept',coefficient(1),'fitSlope',coefficient(2), ...
    'centralFitVelocity',fit,'fitRemainder',u-fit, ...
    'gasWallVelocity',gasU','actualWallSlip',gasU'-cfg.wallVelocity, ...
    'extrapolatedWallSlip',wallFit'-cfg.wallVelocity, ...
    'gasMinusExtrapolatedWallVelocity',(gasU-wallFit)');
d.centralFitRelativeL2 = norm(u(central)-fit(central))/max(norm(u(central)),1e-30);
d.shearStressRelativeSpan = d.shearStressSpan/max(abs(d.meanShearStress),1e-30);
f = figure('Name',sprintf('Couette physics checks, delta=%g',cfg.delta), ...
    'Color','w','Position',[60 60 1400 800]);
if isprop(f,'Theme'); f.Theme = 'light'; end
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile;
plot([wallX(1);x;wallX(2)],[gasU(1);u;gasU(2)],'LineWidth',1.4); hold on;
plot(wallX,wallFit,'--','LineWidth',1.1);
plot(wallX,cfg.wallVelocity,'ks','MarkerFaceColor','k');
labels = {'DG gas velocity','central linear fit','wall velocity'};
if ~isempty(r.reference)
    plot(x,r.reference.cellVelocity,':','LineWidth',1.2);
    labels{end+1} = 'BGK reference (cell average)';
end
xlabel('x'); ylabel('u_y'); title('Tangential velocity and wall slip');
legend(labels,'Location','best'); grid on;

nexttile;
plot(x,q(:,2),'LineWidth',1.2); yline(0,'k:');
span = max(1e-12,1.2*max(abs(q(:,2)))); ylim([-span span]);
xlabel('x'); ylabel('u_x'); title(sprintf('No penetration: max |u_x| = %.2e',d.maxNormalVelocity)); grid on;

nexttile;
plot(x,r.shearStress,'LineWidth',1.3); hold on;
plot(r.mesh.edges,r.diagnostics.faceFlux(:,2),'--','LineWidth',1.1);
span = max([0.05*abs(d.meanShearStress),d.shearStressSpan,1e-8]);
ylim(d.meanShearStress+[-span span]);
xlabel('x'); ylabel('tau_{xy} = <v_x v_y H>_{M_0}');
title(sprintf('Constant shear: relative span = %.2e',d.shearStressRelativeSpan));
legend('cell moment','numerical face flux','Location','best'); grid on;

nexttile;
plot(x,q(:,[1 4]),'LineWidth',1.2);
span = max(1e-12,1.2*max(abs(q(:,[1 4])),[],'all')); ylim([-span span]);
xlabel('x'); ylabel('perturbation'); title('rho = 1 + rho'', T = 1 + T''');
legend('rho''','T''','Location','best'); grid on;

nexttile;
plot(x,d.fitRemainder,'LineWidth',1.3); hold on;
labels = {'DG minus central fit'};
if ~isempty(r.reference)
    % Use the SAME fitted line for DG and reference, to expose wall curvature.
    plot(x,r.reference.cellVelocity-fit,'--','LineWidth',1.1);
    labels{end+1} = 'BGK reference minus central fit';
end
yline(0,'k:');
xline(d.centralFitInterval(1),':','HandleVisibility','off');
xline(d.centralFitInterval(2),':','HandleVisibility','off');
xlabel('x'); ylabel('u_y - u_{central fit}');
title('Departure from central fit');
legend(labels,'Location','best'); grid on;

nexttile;
% Show the first 15% of the gap, including the GAS trace exactly at the wall.
distance = (x-wallX(1))/cfg.gap; near = distance<=0.15;
xx = [0;distance(near)]; uu = [gasU(1);u(near)];
plot(xx,uu,'o-','MarkerSize',3,'LineWidth',1.2); hold on;
plot([0 0.15],coefficient(1)+coefficient(2)*(wallX(1)+[0 0.15]*cfg.gap),'--','LineWidth',1.1);
plot(0,cfg.wallVelocity(1),'ks','MarkerFaceColor','k');
plot([0 0],[cfg.wallVelocity(1) gasU(1)],'k:','HandleVisibility','off');
xlim([-0.005 0.15]); xlabel('(x - x_{left wall}) / gap'); ylabel('u_y');
title(sprintf('Left wall: actual gas-wall slip = %.3e',d.actualWallSlip(1)));
legend('gas velocity and wall trace','central-fit extrapolation','wall velocity','Location','best'); grid on;
sgtitle(sprintf('Couette: x-normal, delta=%g, N=%d; central fit is empirical, not an exact layer decomposition', ...
    cfg.delta,cfg.cells));
drawnow;

fprintf('Couette delta=%g: max |ux| %.3e, max |rho''| %.3e, max |T''| %.3e\n', ...
    cfg.delta,d.maxNormalVelocity,d.maxDensityPerturbation,d.maxTemperaturePerturbation);
fprintf('Shear mean %.9g, relative span %.3e; central-fit relative L2 %.3e\n', ...
    d.meanShearStress,d.shearStressRelativeSpan,d.centralFitRelativeL2);
fprintf('Actual gas-wall slip [left right]: %.9g %.9g\n',d.actualWallSlip);
fprintf('Central-fit extrapolated slip [left right]: %.9g %.9g\n',d.extrapolatedWallSlip);
if saveOutput
    folder = fullfile(cfg.outputDir,'couette');
    if ~exist(folder,'dir'); mkdir(folder); end
    exportgraphics(f,fullfile(folder,'couette_physics.png'),'Resolution',150);
    writematrix([x,u,q(:,2),r.shearStress,q(:,1),q(:,4),fit,d.fitRemainder], ...
        fullfile(folder,'couette_physics.dat'),'Delimiter','tab');
    save(fullfile(folder,'couette_physics.mat'),'d');
    fid = fopen(fullfile(folder,'couette_physics.txt'),'w'); assert(fid>=0,'Cannot write Couette checks.');
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid,'delta=%.9g\ncells=%d\n',cfg.delta,cfg.cells);
    fprintf(fid,'max_abs_ux=%.12e\nmax_abs_rho_pert=%.12e\nmax_abs_T_pert=%.12e\n', ...
        d.maxNormalVelocity,d.maxDensityPerturbation,d.maxTemperaturePerturbation);
    fprintf(fid,'shear_mean=%.12e\nshear_relative_span=%.12e\ncentral_fit_relative_L2=%.12e\n', ...
        d.meanShearStress,d.shearStressRelativeSpan,d.centralFitRelativeL2);
    fprintf(fid,'actual_gas_wall_slip=%.12e %.12e\n',d.actualWallSlip);
    fprintf(fid,'central_fit_extrapolated_slip=%.12e %.12e\n',d.extrapolatedWallSlip);
    fprintf(fid,'fit_interval=%.9g %.9g\nfit_slope=%.12e\n',d.centralFitInterval,d.fitSlope);
    fprintf(fid,'The central fit is empirical; its remainder is not an exact half-space Knudsen-layer function.\n');
    fprintf(fid,'columns: x uy ux tau_xy rho_pert T_pert uy_central_fit uy_minus_fit\n');
end
end
