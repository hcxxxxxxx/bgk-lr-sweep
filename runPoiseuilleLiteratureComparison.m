function comparison = runPoiseuilleLiteratureComparison
%RUNPOISEUILLELITERATURECOMPARISON Compare 1D2V BGK Qp with published data.
% Sharipov, J. Micromech. Microeng. 9 (1999), 394-401, Table 1, p.398.
% Use ONLY a/b=0 (infinite parallel plates), not the rectangular/square data.
% Table columns: present work = S-model; work [36] = BGK; work [29] =
% Boltzmann. These published values are distinct from our integral reference.
%
% Paper: c0=sqrt(2*k*T/m), full height a, delta=a*p/(mu*c0),
% lambda=(mu/p)*sqrt(pi*k*T/(2*m)), Kn=lambda/a=sqrt(pi)/(2*delta).
% Paper driving: upsilon=(a/p)*dp/dx, Q=-Qp*upsilon, Q=2*mean(u).
% Code source: 2*A*vy. With gap=1, upsilon=-2*A, so
% Qp=integral(u_y dx)/A. A=1 is a unit linear response, not finite-Mach flow.
% In the general code geometry Qp=flowRate/(A*gap^2) and
% delta_paper=delta_code*gap; here gap=1 throughout.
%
% Reproduce with: comparison = runPoiseuilleLiteratureComparison;
% Existing saved channel cases and solver defaults are not overwritten.
root = fileparts(mfilename('fullpath'));
folder = fullfile(root,'results','channel_1d2v','literature_poiseuille');
if ~exist(folder,'dir'); mkdir(folder); end
diary(fullfile(folder,'run.log')); diaryCleanup = onCleanup(@() diary('off'));

% Manually transcribed from Table 1; checked against the rendered PDF page.
delta = [0.1; 0.5; 1; 2; 5; 10];
paperBGK = [2.0327; 1.6022; 1.5392; 1.5957; 1.9922; 2.7673];
paperSModel = [2.038; 1.613; 1.549; 1.606; 2.000; 2.772];
paperBoltzmann = [1.9301; 1.5608; 1.5090; 1.5666; 1.9546; 2.7291];
count = numel(delta);
coarseQp = zeros(count,1); fineQp = coarseQp;
coarseIterations = coarseQp; fineIterations = coarseQp;
fineResidual = coarseQp;
cfg = struct('flow','Poiseuille','gap',1,'acceleration',1, ...
    'wallVelocity',[0 0],'initial',0,'tol',1e-10,'stopCriterion','equation','maxIter',10000, ...
    'normalQuadrature','halfLegendre','checkEvery',10,'reportEvery',1000, ...
    'reference',false,'plot',false,'save',false);
for k = 1:count
    cfg.delta = delta(k);
    cfg.cells = 128; cfg.velocityOrder = [256 8];
    coarse = mono1D2V(cfg);
    cfg.cells = 256; cfg.velocityOrder = [512 8];
    fine = mono1D2V(cfg);
    assert(coarse.converged && fine.converged,'Solver did not converge.');
    assert(max(abs(fine.nodalMoments(:,[1 2 4])),[],'all')<1e-9, ...
        'Unexpected normal velocity, density or temperature perturbation.');
    coarseQp(k) = coarse.flowRate/(cfg.acceleration*cfg.gap^2);
    fineQp(k) = fine.flowRate/(cfg.acceleration*cfg.gap^2);
    coarseIterations(k) = coarse.iterations;
    fineIterations(k) = fine.iterations;
    fineResidual(k) = fine.diagnostics.kinetic;
end
refinementPercent = 100*abs(fineQp-coarseQp)./abs(fineQp);
assert(all(refinementPercent<0.01),'Qp refinement change exceeds 0.01 percent.');
bgkDifferencePercent = 100*(fineQp-paperBGK)./paperBGK;
sModelDifferencePercent = 100*(fineQp-paperSModel)./paperSModel;
boltzmannDifferencePercent = 100*(fineQp-paperBoltzmann)./paperBoltzmann;
Kn = sqrt(pi)./(2*delta);
comparison = table(delta,Kn,coarseQp,fineQp,paperBGK,paperSModel,paperBoltzmann, ...
    bgkDifferencePercent,sModelDifferencePercent,boltzmannDifferencePercent, ...
    refinementPercent,coarseIterations,fineIterations,fineResidual);
writetable(comparison,fullfile(folder,'comparison.csv'));

% Also convert the original saved a=0.01 examples rather than replacing them.
oldFiles = {fullfile(root,'results','channel_1d2v','poiseuille','channel.mat'); ...
    fullfile(root,'results','channel_1d2v','delta10','poiseuille','channel.mat')};
savedDelta=[]; savedAcceleration=[]; savedFlow=[]; savedQp=[];
savedBGK=[]; savedBGKDifferencePercent=[];
for k = 1:numel(oldFiles)
    if ~isfile(oldFiles{k}); continue; end
    data = load(oldFiles{k},'r'); r = data.r;
    row = find(delta==r.config.delta*r.config.gap,1);
    qp = r.flowRate/(r.config.acceleration*r.config.gap^2);
    savedDelta(end+1,1) = delta(row); %#ok<AGROW>
    savedAcceleration(end+1,1) = r.config.acceleration; %#ok<AGROW>
    savedFlow(end+1,1) = r.flowRate; %#ok<AGROW>
    savedQp(end+1,1) = qp; %#ok<AGROW>
    savedBGK(end+1,1) = paperBGK(row); %#ok<AGROW>
    savedBGKDifferencePercent(end+1,1) = 100*(qp-paperBGK(row))/paperBGK(row); %#ok<AGROW>
end
savedComparison = table(savedDelta,savedAcceleration,savedFlow,savedQp, ...
    savedBGK,savedBGKDifferencePercent);
writetable(savedComparison,fullfile(folder,'original_cases.csv'));
provenance = struct('title','Non-isothermal gas flow through rectangular microchannels', ...
    'author','Felix Sharipov','year',1999,'journal','J. Micromech. Microeng. 9, 394-401', ...
    'table','Table 1, printed page 398, a/b=0 columns', ...
    'paperBGKSource','Work [36]: S. K. Loyalka, Phys. Fluids 17 (1974), 1053', ...
    'paperSModelSource','Present-work S-model column', ...
    'paperBoltzmannSource','Work [29]: Loyalka and Hickey, J. Vac. Sci. Technol. A 9 (1991), 158', ...
    'normalization','Qp=flowRate/(acceleration*gap^2); gap=1; source=2*acceleration*vy', ...
    'geometry','Fully diffuse stationary isothermal infinite parallel plates; a/b=0', ...
    'note','Paper entries have finite printed precision; differences to other collision models are not BGK discretization errors.');
save(fullfile(folder,'comparison.mat'),'comparison','savedComparison','provenance','cfg');

fig = figure('Name','Poiseuille published Qp comparison','Color','w', ...
    'Position',[100 100 1350 440]);
if isprop(fig,'Theme'); fig.Theme='light'; end
tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
semilogx(delta,fineQp,'o-','LineWidth',1.4);
semilogx(delta,paperBGK,'ks','LineWidth',1.1,'MarkerSize',8);
semilogx(delta,paperSModel,'^--','LineWidth',1.1);
semilogx(delta,paperBoltzmann,'d:','LineWidth',1.1);
set(gca,'XScale','log'); grid on; xlabel('delta'); ylabel('Q_p');
title('Infinite parallel plates: a/b = 0');
legend('Our BGK DG','Table 1: BGK [36]','Table 1: S-model','Table 1: Boltzmann [29]', ...
    'Location','best');
nexttile; hold on;
semilogx(delta,bgkDifferencePercent,'o-','LineWidth',1.4);
semilogx(delta,sModelDifferencePercent,'^--','LineWidth',1.1);
semilogx(delta,boltzmannDifferencePercent,'d:','LineWidth',1.1);
yline(0,'k:'); set(gca,'XScale','log'); grid on;
xlabel('delta'); ylabel('100 (our Q_p / table Q_p - 1) [%]');
title('Differences, separated by collision model');
legend('vs BGK [36]','vs S-model','vs Boltzmann [29]','Location','best');
nexttile;
loglog(delta,refinementPercent,'o-','LineWidth',1.4); grid on;
xlabel('delta'); ylabel('Q_p refinement change [%]');
title({'Resolution check','128 / 256x8 to 256 / 512x8'});
sgtitle('Sharipov (1999), Table 1, p.398: Poiseuille flow coefficient');
exportgraphics(fig,fullfile(folder,'comparison.png'),'Resolution',160);

fid = fopen(fullfile(folder,'comparison.txt'),'w'); assert(fid>=0);
fileCleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s, %s (%d), %s.\n%s.\n',provenance.author,provenance.title, ...
    provenance.year,provenance.journal,provenance.table);
fprintf(fid,'Paper full height a = code gap = 1. Both use c0=sqrt(2*k*T/m).\n');
fprintf(fid,'Paper delta=a*p/(mu*c0); Kn=lambda/a=sqrt(pi)/(2*delta).\n');
fprintf(fid,'Paper upsilon=(a/p)*dp/dx; source=-upsilon*v_stream.\n');
fprintf(fid,'Code source=2*A*v_stream; hence upsilon=-2*A when gap=1.\n');
fprintf(fid,'Paper Q=2*integral(u dx), Q=-Qp*upsilon; hence Qp=flowRate/A.\n');
fprintf(fid,'Unit linear response A=1; stationary fully diffuse isothermal walls.\n');
fprintf(fid,'P1-DG coarse: Nx=128, Vx=256, Vy=8. Fine: Nx=256, Vx=512, Vy=8.\n');
fprintf(fid,'All runs start from zero, absolute tolerance=1e-10, maxIter=10000.\n');
fprintf(fid,'Positive/negative half Gauss-Legendre normal speeds on [0,8]/[-8,0]; tangential Gauss-Hermite.\n');
fprintf(fid,'Published BGK values are transcribed, not our own characteristic-integral reference.\n');
fprintf(fid,'BGK source: %s.\nS-model source: %s.\nBoltzmann source: %s.\n', ...
    provenance.paperBGKSource,provenance.paperSModelSource,provenance.paperBoltzmannSource);
fprintf(fid,'%s\n',provenance.note);
fprintf(fid,'\ndelta our_Qp paper_BGK paper_Smodel paper_Boltzmann BGK_difference_percent refinement_percent\n');
fprintf(fid,'%g %.10f %.4f %.4f %.4f %+.8f %.8f\n', ...
    [delta fineQp paperBGK paperSModel paperBoltzmann bgkDifferencePercent refinementPercent]');
fprintf(fid,'\nOriginal saved a=0.01 cases:\n');
fprintf(fid,'delta acceleration flowRate Qp paper_BGK BGK_difference_percent\n');
fprintf(fid,'%g %g %.12g %.10f %.4f %+.8f\n',table2array(savedComparison)');
disp(comparison);
fprintf('Comparison saved to %s\n',folder);
end
