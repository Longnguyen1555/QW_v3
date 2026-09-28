function cfg = default_config()
%DEFAULT_CONFIG Central configuration for the complete QW-MOAP workflow.

    c = physical_constants();

    %% Run control
    cfg.run.task = 'single';
    cfg.run.verbose = true;
    cfg.run.debug_moap_direct = false;

    %% Material: GaAs
    cfg.material.name = 'GaAs';
    cfg.material.mstar_rel = 0.067;
    cfg.material.eps_static = 13.18;
    cfg.material.eps_high = 10.89;
    cfg.material.LO_phonon_meV = 36.25;
    cfg.material.sound_speed_mps = 5.22e3;
    cfg.material.mass_density_kgm3 = 5317;
    cfg.material.piezo_kappa2 = 0.006;
    cfg.material.e14_Cm2 = 0.16;

    cfg.material.eps_poisson = 13.18;

    %% Quantum-well confinement
    cfg.structure.potential_model = 'manning_sech';

    cfg.structure.U0_meV = 220.0;
    cfg.structure.Lz_nm = 12.0;
    cfg.structure.manning_prefactor = 6.0;

    cfg.structure.alpha = 0.30;
    cfg.structure.domain_nm = 60.0;
    cfg.structure.auto_domain = false;
    cfg.structure.domain_factor = 3.0;
    cfg.structure.Nz = 1502;
    cfg.structure.n_states = 6;

    %% Delta doping
    cfg.doping.Nd_sheet_cm2 = 1.0e13;
    cfg.doping.width_nm = 2.0;

    %% Static fields
    cfg.fields.B_T = 0.0;
    cfg.fields.E_kVcm = 0.0;
    cfg.fields.kx_inv_m = 0.0;

    %% Temperature
    cfg.temperature_K = 90.0;

    %% Schrodinger-Poisson controls
    cfg.sp.max_outer = 100;
    cfg.sp.max_inner = 500;
    cfg.sp.mix = 0.40; % stable PDF fixed-point mixing for configured subbands
    cfg.sp.tol_VH_meV = 1.0e-3;
    cfg.sp.tol_VH_relative = 1.0e-3; % PDF relative U_H iteration criterion
    cfg.sp.tol_EF_meV = 1.0e-3;
    cfg.sp.poisson_boundary = 'dirichlet_zero';

    %% Laser
    cfg.laser.E0_kVcm = 4.5;
    cfg.laser.a0_mode = 'constant'; %'dynamic'
    cfg.laser.a0_nm = 7.5;

    %% OAP model
    cfg.oap.model = 'direct_q_integral';
    cfg.oap.transitions = [1 2];
    cfg.oap.photon_orders = [1 2 3];
    cfg.oap.mechanisms = {'optical', 'piezoelectric'};
    cfg.oap.photon_energy_meV = 7.5:0.01:100.0;
    cfg.oap.population_mode = 'maxwell_boltzmann';
    cfg.oap.include_pauli_blocking = false;

    % q integration: q is in nm^-1 at the user level
    cfg.oap.qz_max_inv_nm = 2.5;
    cfg.oap.qperp_max_inv_nm = 1.5;
    cfg.oap.Nqz = 90;
    cfg.oap.Nqperp = 90;

    cfg.oap.inplane_factor = 'none';
    cfg.oap.landau_initial = 0;
    cfg.oap.landau_final = 0;
    cfg.oap.thermal_cutoff_factor = 1.0;

    cfg.oap.gamma_optical_meV = 0.80;
    cfg.oap.gamma_piezo_meV = 0.60;

    cfg.oap.screening_density_mode = 'fixed';
    cfg.oap.fixed_ne_cm3 = 1.0e18;
    cfg.oap.active_thickness_nm = 10.0;

    cfg.oap.plot_normalization = 'none';
    % Keep the complete configured photon-energy interval on spectrum plots.
    % This is required when a resonance lies on the negative-energy branch.
    cfg.oap.plot_full_energy_range = true;
    cfg.oap.direct.RelTol = 1.0e-9;
    cfg.oap.direct.AbsTol = 0.0;

    % Largest verified truncations below 200.  s_max and v_max remain finite
    % at 199 in the Refined checks; Nhan overflows at eta=18 for a valid
    % D>0 fixture, so eta_max=17 is the largest finite verified value.
    cfg.oap.series.s_max = 199;   % screening index s = 0:s_max (van), 1:s_max (Nhan)
    cfg.oap.series.eta_max = 17;  % HƯỚNG DẪN screening-series index eta
    cfg.oap.series.v_max = 199;   % Refined PDF Fermi index v = 1:v_max
    cfg.oap.series.evaluation = 'auto'; % series = checked Eq. (32), quadrature = Eq. (28)
    cfg.oap.series.enforce_fermi_condition = true;
    cfg.oap.normalization_area_m2 = 1.0;

    %% Sweep
    cfg.sweep.parameter = 'B_T';
    cfg.sweep.values = [5 10 15];

    %% Output
    cfg.output.directory = fullfile(pwd, 'results');
    cfg.output.save_figures = true;
    cfg.output.save_csv = true;
    cfg.output.save_mat = true;
    cfg.output.visible = 'on';
    cfg.output.figure_format = 'png';
    cfg.output.figure_dpi = 200;

    %% Constants kept for traceability
    cfg.constants = c;
end
