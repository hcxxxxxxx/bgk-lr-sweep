function result = monoGSIS1D2V(options)
%MONOGSIS1D2V Full-state 1D2V BGK/HDG coupling with kinetic-flux walls.
% r=monoGSIS1D2V; runs Couette and Poiseuille, with S=[1 delta delta delta].
% r=monoGSIS1D2V(struct('flow','Poiseuille','delta',10,'plot',false));
% Qsrc=Qmacro is returned directly, without an H lift or small mixing factor.
% HOT derivatives include actual upwind face moments. Wall traces follow
% F_HDG=F_kinetic+J*(U-Ukinetic); the mean density gauge is fixed to zero.
% A seven-moment weak load correction preserves the kinetic-DG fixed point.
% options.macroParam.kineticConsistency=false disables that extra load only.
% Default stopping matches reference cell-macro relative updates, tol=1e-5.
% diagnostics.kinetic/coupled separately report both equation residuals;
% stopCriterion='equation' optionally restores residual-based stopping.
if nargin<1; options=struct; end
options.iterationMethod='GSIS';
if ~isfield(options,'outputDir')
    options.outputDir=fullfile(fileparts(mfilename('fullpath')),'results','channel_1d2v_gsis');
end
result=mono1D2V(options);
end
