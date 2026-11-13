function plot_su2019_bgk_couette
% Relative to independent values printed in Table 1; all plotted errors are %.
root=fileparts(fileparts(mfilename('fullpath')));
folder=fullfile(root,'results','su2019_bgk_couette');
b=readtable(fullfile(folder,'baseline_summary.csv'));
t=readtable(fullfile(folder,'refinement_summary.csv'));
groups={b(strcmp(b.method,'CIS'),:),b(strcmp(b.method,'GSIS'),:), ...
    t(strcmp(t.method,'CIS'),:),t(strcmp(t.method,'GSIS') & t.cells==128,:), ...
    t(strcmp(t.method,'GSIS') & t.cells==512,:)};
labels={'CIS N128 V32x32 tol1e-5','GSIS N128 V32x32 tol1e-5', ...
    'CIS N128 V32x32 tol1e-8','GSIS N128 V32x32 tol1e-8','GSIS N512 V64x32 tol1e-8'};
styles={'o-','s-','o:','s:','^--'};colors=[0 .35 .75;.8 .2 .1;0 .35 .75;.8 .2 .1;.15 .15 .15];
fields={'wallErrorPercent','gradientErrorPercent','shearErrorPercent'};
titles={'Gas velocity at the wall','Velocity gradient at channel center','Shear stress'};
f=figure('Visible','off','Color','w','Position',[60 60 1400 470]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for metric=1:3
    nexttile;hold on;
    for group=1:numel(groups)
        q=groups{group};loglog(q.delta,max(q.(fields{metric}),1e-12),styles{group}, ...
            'Color',colors(group,:),'LineWidth',1.3,'MarkerSize',5);
    end
    set(gca,'XScale','log','YScale','log');grid on;xlabel('Collision strength delta');
    ylabel('Absolute relative error (%)');title(titles{metric});
    if metric==1;legend(labels,'Location','best','FontSize',8,'Interpreter','none');end
end
sgtitle('BGK Couette vs Su et al. JCP 378 (2019), Table 1; measured from kinetic H');
exportgraphics(f,fullfile(folder,'paper_comparison.png'),'Resolution',180);close(f);
f=figure('Visible','off','Color','w','Position',[60 60 1050 450]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for metric=1:2
    nexttile;hold on;
    for group=1:4
        q=groups{group};y=q.iterations;if metric==2;y=q.relativeOriginalResidual;end
        loglog(q.delta,y,styles{group},'Color',colors(group,:),'LineWidth',1.3,'MarkerSize',5);
    end
    set(gca,'XScale','log','YScale','log');grid on;xlabel('Collision strength delta');
    if metric==1;ylabel('Iterations');title('Cell-macro stopping');legend(labels(1:4),'Location','best','FontSize',8,'Interpreter','none');
    else;ylabel('Original kinetic residual / first-sweep residual');title('Original equation, independent of stopping');end
end
exportgraphics(f,fullfile(folder,'stopping_and_residual.png'),'Resolution',180);close(f);
end
