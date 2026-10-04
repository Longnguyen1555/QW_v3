
clear; clc; close all;

project_root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(project_root, 'config')));
addpath(genpath(fullfile(project_root, 'src')));

cfg = default_config();

%% ========================================================================
% USER SETTINGS
% ========================================================================
% Available tasks:
%   'single'           : one electronic structure + one MOAP spectrum
%   'sweep_spectra'    : overlay spectra while varying one parameter
%   'linewidth_sweep'  : FWHM/HWHM versus one parameter
%   'energy_gap_sweep' : E2-E1 from the self-consistent SP solver
%   'all_demo'         : single + sweep spectra + linewidth sweep
cfg.run.task = 'single';

% Numerical profile: 'quick', 'standard', or 'high_accuracy'
cfg = apply_numerical_profile(cfg, 'standard');

% Main physical parameters
cfg.structure.Lz_nm       = 10.0;
cfg.structure.U0_meV      = 228.0;
cfg.doping.Nd_sheet_cm2   = 1.0e13;
cfg.fields.B_T            = 5.0;
cfg.fields.E_kVcm         = 25.0;
cfg.temperature_K         = 10.0;

% Transition indices are MATLAB 1-based subband indices.
cfg.oap.model             = 'direct_approx'; %lorentzian, van, nhan, direct_q_integral
cfg.oap.transitions       = [1 2];
cfg.oap.photon_orders     = [1 2];
cfg.oap.mechanisms        = {'optical', 'piezoelectric'};%'optical', 'piezoelectric'

% Photon-energy scan
cfg.oap.photon_energy_meV = 20.0:0.01:300.0;

% Sweep used by 'sweep_spectra', 'linewidth_sweep', and 'energy_gap_sweep'.
% For the energy-gap task, B_T may be changed to T_K, Lz_nm, E_kVcm,
% or Nd_sheet_cm2 without changing the task implementation.
% Supported names: 'B_T', 'T_K', 'Lz_nm', 'E_kVcm',
%                  'Nd_sheet_cm2'
cfg.sweep.parameter       = 'T_K';
cfg.sweep.values          = [10 20 30];

% Output
cfg.output.directory      = fullfile(project_root, 'results_approx');%results_van, results_nhan, results_lorentzian
cfg.output.save_figures   = true;
cfg.output.save_csv       = true;
cfg.output.save_mat       = true;
cfg.output.visible        = 'on';

%% Run
results = run_qw_moap(cfg);

fprintf('\n[MOAP] Completed successfully.\n');
fprintf('[MOAP] Results folder: %s\n', cfg.output.directory);
