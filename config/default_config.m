function cfg = default_config()
%DEFAULT_CONFIG Central configuration for the complete QW-MOAP workflow.

    c = physical_constants();

    %% Run control
    cfg.run.task = 'single';
    cfg.run.verbose = true;

    %% Material: GaAs
    cfg.material.name = 'GaAs';
    cfg.material.mstar_rel = 0.067;
    cfg.material.eps_static = 13.18;
    cfg.material.eps_high = 10.89;
    cfg.material.LO_phonon_meV = 36.25;
    cfg.material.sound_speed_mps = 5.22e3;
    cfg.material.mass_density_kgm3 = 5317;
    cfg.material.piezo_kappa2 = 0.006;
    cfg.material.e14_Cm2 = 0.16;  % optional tensor-based estimate

    %% Anharmonic confinement:
    % U_A(z) = U0 * (z/Lz)^2 * [alpha*(z/Lz)^6 - 1]
    cfg.structure.U0_meV = 228.0;
    cfg.structure.Lz_nm = 5.0;
    cfg.structure.alpha = 0.30;
    cfg.structure.domain_nm = 15.0;
    cfg.structure.auto_domain = true;
    cfg.structure.domain_factor = 3.0;
    cfg.structure.Nz = 1200;
    cfg.structure.n_states = 4;

    %% Delta doping
    cfg.doping.Nd_sheet_cm2 = 1.0e13;
    cfg.doping.width_nm = 2.0;

    %% Static fields
    cfg.fields.B_T = 10.0;
    cfg.fields.E_kVcm = 0.0;
    cfg.fields.kx_inv_m = 0.0;

    %% Temperature
    cfg.temperature_K = 300.0;

    %% Schrodinger-Poisson controls
    cfg.sp.max_iter = 300;
    cfg.sp.mix = 0.20;
    cfg.sp.tol_potential_rel = 1.0e-4;
    cfg.sp.tol_EF_meV = 1.0e-6;
    cfg.sp.minimum_scale_meV = 1.0e-3;
    cfg.sp.poisson_boundary = 'dirichlet_zero';

    %% Laser
    cfg.laser.E0_kVcm = 4.5;
    cfg.laser.a0_mode = 'constant';  % 'constant' or 'dynamic'
    cfg.laser.a0_nm = 7.5;

    %% OAP model
    % 'direct_q_integral' evaluates the transition-probability integral
    % before the Taylor-Maclaurin expansion. This is the recommended mode.
    % 'analytical_series' implements the truncated Eq. (5)/(7) series.
    cfg.oap.model = 'direct_q_integral';
    cfg.oap.transitions = [1 2];
    cfg.oap.photon_orders = [1 2 3];
    cfg.oap.mechanisms = {'optical', 'piezoelectric'};
    cfg.oap.photon_energy_meV = 2.0:0.10:160.0;
    cfg.oap.population_mode = 'maxwell_boltzmann';
    cfg.oap.include_pauli_blocking = false;

    % q integration: q is in nm^-1 at the user level
    cfg.oap.qz_max_inv_nm = 2.5;
    cfg.oap.qperp_max_inv_nm = 1.5;
    cfg.oap.Nqz = 90;
    cfg.oap.Nqperp = 90;

    % In-plane factor. 'landau' uses Eq. (12), reducing to exp(-u) for N=0.
    % At B=0 the code automatically uses a thermal recoil cutoff.
    cfg.oap.inplane_factor = 'landau';
    cfg.oap.landau_initial = 0;
    cfg.oap.landau_final = 0;
    cfg.oap.thermal_cutoff_factor = 1.0;

    % Collision broadening used to replace the Dirac delta by a Lorentzian.
    cfg.oap.gamma_optical_meV = 0.80;
    cfg.oap.gamma_piezo_meV = 0.60;

    % Debye screening density. 'from_sheet' converts sheet density using
    % the active electronic thickness; 'fixed' uses fixed_ne_cm3.
    cfg.oap.screening_density_mode = 'fixed';
    cfg.oap.fixed_ne_cm3 = 1.0e18;
    cfg.oap.active_thickness_nm = 10.0;

    % Plot normalization does not alter saved raw spectra.
    cfg.oap.plot_normalization = 'global_max'; % 'none' or 'global_max'

    % Truncated analytical-series controls
    cfg.oap.series.s_max = 8;
    cfg.oap.series.eta_max = 5;
    cfg.oap.series.piezo_effective_q_mode = 'debye';

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
