% Overlay MOAP spectra for several magnetic fields.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'config')));
addpath(genpath(fullfile(root,'src')));

cfg = apply_numerical_profile(default_config(),'standard');
cfg.run.task = 'sweep_spectra';
cfg.sweep.parameter = 'B_T';
cfg.sweep.values = [5 10 15];
cfg.output.directory = fullfile(root,'results','spectrum_vs_B');
run_qw_moap(cfg);
