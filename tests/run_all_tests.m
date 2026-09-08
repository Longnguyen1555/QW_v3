function run_all_tests()
%RUN_ALL_TESTS Toolbox-free unit and integration tests.

    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(fullfile(root,'config')));
    addpath(genpath(fullfile(root,'src')));

    fprintf('Running QW-MOAP tests...\n');

    test_wavefunction_normalization();
    test_charge_neutrality();
    test_fwhm_synthetic_peak();
    test_fwhm_boundary_peak();
    test_fwhm_multiple_peaks();
    test_legacy_broadening_width();
    test_fwhm_does_not_change_sp();
    test_resonance_positions();
    test_physics_configuration();
    test_full_pipeline();

    fprintf('ALL TESTS PASSED\n');
end
