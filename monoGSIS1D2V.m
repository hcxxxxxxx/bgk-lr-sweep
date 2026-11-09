function result = monoGSIS1D2V(options)
%MONOGSIS1D2V Reference-style 1D2V BGK/HDG coupled iteration, without H lift.
% r=monoGSIS1D2V; runs Couette and Poiseuille, with S=[1 delta delta delta].
% r=monoGSIS1D2V(struct('flow','Poiseuille','delta',10,'plot',false));
% Only source moments are mixed: Qsrc=Qkin+omega*(Qmacro-Qkin), where
% omega=min(1,1/Kn_cell), Kn_cell=sqrt(pi)/(2*delta*dx).
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
