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
    cfg.material.e14_Cm2 = 0.16;

    cfg.material.eps_poisson = 13.18;

    %% Quantum-well confinement
    cfg.structure.potential_model = 'manning_sech';

    cfg.structure.U0_meV = 220.0;
    cfg.structure.Lz_nm = 12.0;
    cfg.structure.manning_prefactor = 6.0;

    % alpha belonged to the retired polynomial potential and has no
    % source-supported meaning for the Manning/sech model.
    cfg.structure.alpha = NaN;
    cfg.structure.domain_nm = 60.0;
    cfg.structure.auto_domain = false;
    cfg.structure.domain_factor = 3.0;
    cfg.structure.Nz = 1502;
    % Historical SP baseline solved 12 states; keep that basis explicit.
    cfg.structure.n_states = 12;

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
    cfg.sp.mix = 0.50;
    cfg.sp.tol_VH_meV = 1.0e-3;
    cfg.sp.tol_EF_meV = 1.0e-3;
    cfg.sp.poisson_boundary = 'dirichlet_zero';
    cfg.sp.occupation_warning_fraction = 1.0e-3;

    %% Laser
    cfg.laser.E0_kVcm = 4.5;
    cfg.laser.a0_mode = 'dynamic'; %'dynamic'
    cfg.laser.a0_nm = 7.5;

    %% OAP model
    cfg.oap.model = 'direct_q_integral';
    cfg.oap.transitions = [1 2];
    % The supplied detailed derivation supports only 1PA and 2PA.
    cfg.oap.photon_orders = [1 2];
    cfg.oap.mechanisms = {'optical', 'piezoelectric'};
    % Dynamic a0 is undefined at zero photon energy.
    cfg.oap.photon_energy_meV = 2.0:0.010:100.0;
    cfg.oap.population_mode = 'maxwell_boltzmann';
    cfg.oap.include_pauli_blocking = false;

    % q integration: q is in nm^-1 at the user level
    cfg.oap.qz_max_inv_nm = 2.5;
    cfg.oap.qperp_max_inv_nm = 1.5;
    cfg.oap.Nqz = 90;
    cfg.oap.Nqperp = 90;
    % Source derivation: q ~= q_perp, with the qz form factor separated.
    % generalized_3d retains q = sqrt(q_perp^2 + qz^2) for comparison.
    cfg.oap.q_model = 'source_qperp';
    % Controlled quadrature for the momentum left after delta integration.
    cfg.oap.Nkin_perp = 80;
    cfg.oap.fermi_tail_kBT = 40;

    cfg.oap.inplane_factor = 'none';
    cfg.oap.landau_initial = 0;
    cfg.oap.landau_final = 0;
    cfg.oap.thermal_cutoff_factor = 1.0;

    % Optional phenomenological display broadening.  These widths are used
    % only when broadening.mode is 'legacy_lorentzian'; they are not a
    % microscopic linewidth model.
    cfg.oap.broadening.mode = 'none'; % 'none' or 'legacy_lorentzian'
    cfg.oap.gamma_optical_meV = 0.80;
    cfg.oap.gamma_piezo_meV = 0.60;

    cfg.oap.screening_density_mode = 'fixed';
    cfg.oap.fixed_ne_cm3 = 1.0e18;
    cfg.oap.active_thickness_nm = 10.0;

    cfg.oap.plot_normalization = 'global_max';
    cfg.oap.direct.RelTol = 1.0e-9;
    cfg.oap.direct.AbsTol = 0.0;

    cfg.oap.series.s_max = 8;
    cfg.oap.series.eta_max = 5;
    cfg.oap.series.v_max = 20;
    cfg.oap.series.enforce_fermi_condition = true;
    cfg.oap.normalization_area_m2 = 1.0;
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
