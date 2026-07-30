% FWHM/HWHM versus magnetic field.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'config')));
addpath(genpath(fullfile(root,'src')));

cfg = apply_numerical_profile(default_config(),'standard');
cfg.run.task = 'linewidth_sweep';
cfg.sweep.parameter = 'B_T';
cfg.sweep.values = 2:2:20;
cfg.output.directory = fullfile(root,'results','linewidth_vs_B');
run_qw_moap(cfg);
