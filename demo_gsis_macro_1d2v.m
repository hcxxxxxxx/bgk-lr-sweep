function result = demo_gsis_macro_1d2v(options)
%DEMO_GSIS_MACRO_1D2V Run the kinetic-to-synthetic HDG macroscopic interface.
% result=demo_gsis_macro_1d2v;  % Couette and Poiseuille, default 32x32 velocities
% result=demo_gsis_macro_1d2v(struct('delta',10,'velocityOrder',[64 32]));
% One HDG subsolve is performed after each converged kinetic example.
% This demo verifies the macroscopic module; it does not measure acceleration
% of a coupled GSIS outer loop. rho and T are perturbations throughout.
% options are forwarded to mono1D2V (kinetic plot/save are suppressed).
% Results and plots are saved separately in results/gsis_macro_1d2v.
if nargin<1; options=struct; end
options.plot=false; options.save=false;
options.iterationMethod='CIS';
kinetic=mono1D2V(options);
root=fileparts(mfilename('fullpath')); folder=fullfile(root,'results','gsis_macro_1d2v');
if ~exist(folder,'dir'); mkdir(folder); end
if isfield(kinetic,'couette'); names={'couette','poiseuille'};
else; names={lower(kinetic.config.flow)}; end
result=struct;
window=figure('Name','1D2V synthetic HDG macroscopic subsolve', ...
    'Color','w','Position',[100 100 1350 420*numel(names)]);
if isprop(window,'Theme'); window.Theme='light'; end
tiledlayout(window,numel(names),3,'TileSpacing','compact','Padding','compact');
metrics=zeros(numel(names),9);
for i=1:numel(names)
    name=names{i};
    if isfield(kinetic,'couette'); k=kinetic.(name); else; k=kinetic; end
    [macro,aux,hdg,hot]=gsis_macro_from_kinetic_1d2v(k);
    result.(name)=struct('kinetic',k,'macro',macro,'aux',aux,'hdg',hdg,'hot',hot);
    macroUy=macro.cellAverage(:,3); kineticUy=k.cellMoments(:,3);
    relativeChange=norm(macroUy-kineticUy)/norm(kineticUy);
    nodalChange=norm(macro.uy-k.nodalMoments(:,3))/norm(k.nodalMoments(:,3));
    metrics(i,:)=[k.config.delta,k.iterations,relativeChange,hdg.localResidual, ...
        hdg.traceResidual,hdg.fluxContinuityResidual,hdg.conservationResidual, ...
        max(abs(hot.velocityMomentDefect)),nodalChange];
    fprintf('HDG %s: kinetic sweeps=%d; uy difference=%.6e; local/trace/jump/balance=%.3e %.3e %.3e %.3e\n', ...
        name,k.iterations,relativeChange,metrics(i,4:7));
    fprintf('  These are kinetic-CIS sweep counts followed by ONE HDG solve, not GSIS iteration counts.\n');
    x=k.mesh.center;
    nexttile; plot(x,kineticUy,'LineWidth',1.4); hold on;
    plot(x,macroUy,'--','LineWidth',1.4);
    if k.config.reference
        plot(x,k.reference.cellVelocity,':','LineWidth',1.3);
        legendNames={'kinetic DG','one synthetic HDG solve','BGK integral reference'};
    else
        legendNames={'kinetic DG','one synthetic HDG solve'};
    end
    xlabel('x'); ylabel('u_y'); title(sprintf('%s: delta=%g',name,k.config.delta));
    legend(legendNames,'Location','best'); grid on;
    % The conserved discrete quantity is the numerical face flux. The local
    % constitutive stress may jump between DG cells; save it in profile.csv.
    facePxy=[-hdg.outwardFlux(3,1);hdg.outwardFlux(7,:)']/2;
    nexttile; plot(k.mesh.edges,facePxy,'LineWidth',1.3); hold on;
    plot(k.mesh.edges,k.diagnostics.faceFlux(:,2),'--','LineWidth',1.3);
    xlabel('x'); ylabel('shear momentum flux'); grid on;
    legend('HDG numerical flux','kinetic numerical flux','Location','best');
    title('Numerical face flux');
    axisHandle=gca; axisHandle.YAxis.Exponent=0;
    allFlux=[facePxy;k.diagnostics.faceFlux(:,2)];
    if strcmpi(name,'couette') || max(allFlux)-min(allFlux)<1e-3*max(abs(allFlux))
        centerFlux=mean(allFlux);
        padding=max([0.01*abs(centerFlux),0.6*(max(allFlux)-min(allFlux)),1e-12]);
        ylim(centerFlux+[-padding padding]);
    end
    ytickformat('%.7g');
    nexttile;
    plot(x,macroUy-kineticUy,'LineWidth',1.3); grid on;
    xlabel('x'); ylabel('u_y^{HDG} - u_y^{kinetic}');
    title({'Cell-mean velocity change',sprintf('Relative L2 = %.3e',relativeChange)});
    axisHandle=gca; axisHandle.YAxis.Exponent=0; ytickformat('%.1e');
    profile=[hdg.cache.nodes,k.nodalMoments,macro.nodal,aux.Lx,aux.Ly,aux.E, ...
        aux.sigma_xx,aux.sigma_xy,aux.qx,hot.derivative];
    profileNames={'x','kinRho','kinUx','kinUy','kinT','hdgRho','hdgUx','hdgUy','hdgT', ...
        'Lx','Ly','E','sigma_xx','sigma_xy','qx','Pi_xx','Pi_xy','Theta_x'};
    writetable(array2table(profile,'VariableNames',profileNames),fullfile(folder,[name,'_profile.csv']));
    faceTable=array2table([k.mesh.edges,facePxy,k.diagnostics.faceFlux(:,2)], ...
        'VariableNames',{'x','hdgShearMomentumFlux','kineticShearMomentumFlux'});
    writetable(faceTable,fullfile(folder,[name,'_face_flux.csv']));
end
sgtitle('Kinetic DG and one synthetic P1-HDG macroscopic solve');
for axis=findall(window,'Type','axes')'; axis.Toolbar.Visible='off'; end
drawnow;
exportgraphics(window,fullfile(folder,'macro_comparison.png'),'Resolution',150);
summary=array2table(metrics,'VariableNames',{'delta','kineticCISIterations', ...
    'relativeCellUyChange','localResidual','traceResidual','fluxJump','conservationResidual', ...
    'velocityMomentDefect','relativeNodalUyChange'});
summary.flow=names'; summary=movevars(summary,'flow','Before',1);
writetable(summary,fullfile(folder,'summary.csv'));
save(fullfile(folder,'macro_examples.mat'),'result','summary');
fprintf('Macroscopic examples saved to %s\n',folder);
end
