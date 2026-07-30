% Main_quick_test.m
% Small deterministic smoke test. It should complete much faster than Main.m.

clear; clc; close all;

project_root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(project_root, 'config')));
addpath(genpath(fullfile(project_root, 'src')));

cfg = default_config();
cfg = apply_numerical_profile(cfg, 'quick');
cfg.run.task = 'single';
cfg.output.directory = fullfile(project_root, 'results', 'quick_test');
cfg.output.visible = 'off';
cfg.output.save_figures = true;
cfg.output.save_csv = true;
cfg.output.save_mat = true;

results = run_qw_moap(cfg);

assert(results.single.sp.converged, ...
    'Quick test failed: Schrodinger-Poisson solver did not converge.');
assert(all(isfinite(results.single.spectrum.total_raw)), ...
    'Quick test failed: spectrum contains non-finite values.');
assert(max(results.single.spectrum.total_raw) > 0, ...
    'Quick test failed: spectrum is identically zero.');

fprintf('[QUICK TEST] PASS\n');
fprintf('[QUICK TEST] E1..E4 (meV): ');
fprintf('%.6f ', results.single.sp.E_meV);
fprintf('\n');
