% Optical-phonon-only spectrum.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'config')));
addpath(genpath(fullfile(root,'src')));

cfg = apply_numerical_profile(default_config(),'standard');
cfg.run.task = 'single';
cfg.oap.mechanisms = {'optical'};
cfg.output.directory = fullfile(root,'results','optical_only');
run_qw_moap(cfg);
