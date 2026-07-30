% Overlay MOAP spectra for several temperatures.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'config')));
addpath(genpath(fullfile(root,'src')));

cfg = apply_numerical_profile(default_config(),'standard');
cfg.run.task = 'sweep_spectra';
cfg.sweep.parameter = 'T_K';
cfg.sweep.values = [77 150 300];
cfg.output.directory = fullfile(root,'results','spectrum_vs_T');
run_qw_moap(cfg);
