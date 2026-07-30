function run_all_tests()
%RUN_ALL_TESTS Toolbox-free unit and integration tests.

    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(fullfile(root,'config')));
    addpath(genpath(fullfile(root,'src')));

    fprintf('Running QW-MOAP tests...\n');

    test_wavefunction_normalization();
    test_charge_neutrality();
    test_resonance_positions();
    test_full_pipeline();

    fprintf('ALL TESTS PASSED\n');
end
