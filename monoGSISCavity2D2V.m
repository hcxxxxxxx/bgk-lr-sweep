function r=monoGSISCavity2D2V(options)
% Genuine 2D2V cavity GSIS, Q1-HDG macro source return, no distribution lift.
% Both original and coupled kinetic residuals are recorded independently.
% Default algorithm follows GSIS_Gray_20261110: unit stabilization,
% HOT-only internal load, full macro return, cell-integrated macro stopping.
% The physical model remains 2V BGK, rather than the archive's Callaway model.
% Residual matching and relaxation are available only by explicit request.
if nargin<1;options=struct;end
options.iterationMethod='GSIS';
if ~isfield(options,'macroParam');options.macroParam=struct;end
if ~isfield(options.macroParam,'kineticConsistency');options.macroParam.kineticConsistency=false;end
if ~isfield(options.macroParam,'relaxation');options.macroParam.relaxation=1;end
if ~isfield(options,'tol');options.tol=1e-5;end
if ~isfield(options,'outputDir')
    options.outputDir=fullfile(fileparts(mfilename('fullpath')),'results','cavity_gsis_gray_defaults_2d2v');
end
r=monoCavity2D2V(options);
end
