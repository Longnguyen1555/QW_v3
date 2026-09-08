function test_fwhm_multiple_peaks()
%TEST_FWHM_MULTIPLE_PEAKS Select the dominant peak without merging peaks.

    x = -8:0.001:8;
    sigma = 0.35;
    y = exp(-0.5*((x+3)/sigma).^2) + ...
        0.8*exp(-0.5*((x-3)/sigma).^2);
    m = profile_fwhm(x, y);

    expected = 2*sqrt(2*log(2))*sigma;
    assert(abs(m.peak_x + 3) < 2e-3, ...
        'The dominant interior peak was not selected.');
    assert(abs(m.fwhm - expected) < 2e-3, ...
        'FWHM crossed a separate physical peak.');
    fprintf('  PASS: multiple FWHM peaks\n');
end
