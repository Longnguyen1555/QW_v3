% Main.m
% Full MATLAB workflow for nonlinear optical absorption power (MOAP)
% in a Si delta-doped anharmonic GaAs quantum well.
%
% The default run:
%   1) solves the self-consistent Schrodinger-Poisson problem;
%   2) computes optical-phonon and piezoelectric-phonon spectra;
%   3) plots 1PA, 2PA, 3PA contributions;
%   4) saves MAT, CSV and PNG outputs in results/.
%
% Edit only the USER SETTINGS block below for normal use.

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
%   'all_demo'         : single + sweep spectra + linewidth sweep
cfg.run.task = 'all_demo';

% Numerical profile: 'quick', 'standard', or 'high_accuracy'
cfg = apply_numerical_profile(cfg, 'standard');

% Main physical parameters
cfg.structure.Lz_nm       = 5.0;
cfg.structure.alpha       = 0.30;
cfg.structure.U0_meV      = 220.0;
cfg.doping.Nd_sheet_cm2   = 1.0e13;
cfg.doping.width_nm       = 2.0;
cfg.fields.B_T            = 10.0;
cfg.fields.E_kVcm         = 0.0;
cfg.temperature_K         = 90.0;

% Transition indices are MATLAB 1-based subband indices.
cfg.oap.model             = 'direct_q_integral';
cfg.oap.transitions       = [1 2];
cfg.oap.photon_orders     = [1 2 3];
cfg.oap.mechanisms        = {'optical', 'piezoelectric'};

% Photon-energy scan
cfg.oap.photon_energy_meV = 2.0:0.10:160.0;

% Sweep used by 'sweep_spectra' and 'linewidth_sweep'
% Supported names: 'B_T', 'T_K', 'Lz_nm', 'E_kVcm',
%                  'Nd_sheet_cm2', 'alpha', 'U0_meV'
cfg.sweep.parameter       = 'B_T';
cfg.sweep.values          = [5 10 15];

% Output
cfg.output.directory      = fullfile(project_root, 'results');
cfg.output.save_figures   = true;
cfg.output.save_csv       = true;
cfg.output.save_mat       = true;
cfg.output.visible        = 'on';

%% Run
results = run_qw_moap(cfg);

fprintf('\n[MOAP] Completed successfully.\n');
fprintf('[MOAP] Results folder: %s\n', cfg.output.directory);
