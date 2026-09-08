function test_fwhm_boundary_peak()
%TEST_FWHM_BOUNDARY_PEAK Boundary maxima must not receive a fake width.

    x = 0:0.01:10;
    y = exp(-x);
    m = profile_fwhm(x, y);

    assert(m.is_boundary_peak, 'Boundary peak was not flagged.');
    assert(isnan(m.fwhm) && isnan(m.hwhm), ...
        'Boundary peak incorrectly received a finite linewidth.');
    fprintf('  PASS: boundary FWHM peak\n');
end
