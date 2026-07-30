% Overlay MOAP spectra for several characteristic well widths.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'config')));
addpath(genpath(fullfile(root,'src')));

cfg = apply_numerical_profile(default_config(),'standard');
cfg.run.task = 'sweep_spectra';
cfg.sweep.parameter = 'Lz_nm';
cfg.sweep.values = [4.5 5.0 5.5];
cfg.output.directory = fullfile(root,'results','spectrum_vs_Lz');
run_qw_moap(cfg);
