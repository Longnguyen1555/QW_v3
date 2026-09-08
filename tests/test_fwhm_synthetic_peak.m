function test_fwhm_synthetic_peak()
%TEST_FWHM_SYNTHETIC_PEAK Verify crossings without assuming one line shape.

    x = -10:0.001:10;

    gamma = 0.8;
    lorentzian = gamma/pi ./ (x.^2 + gamma^2);
    m_lorentzian = profile_fwhm(x, lorentzian);
    assert(abs(m_lorentzian.fwhm - 2*gamma) < 2e-3, ...
        'Lorentzian FWHM crossing test failed.');

    sigma = 0.7;
    gaussian = exp(-0.5*((x-1.0)/sigma).^2);
    m_gaussian = profile_fwhm(x, gaussian);
    expected_gaussian = 2*sqrt(2*log(2))*sigma;
    assert(abs(m_gaussian.fwhm - expected_gaussian) < 2e-3, ...
        'Gaussian FWHM crossing test failed.');
    fprintf('  PASS: synthetic FWHM peaks\n');
end
