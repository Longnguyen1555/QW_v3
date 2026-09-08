function test_legacy_broadening_width()
%TEST_LEGACY_BROADENING_WIDTH Preserve the historical gamma/order widths.

    x = -20:0.001:20;
    gamma = 0.8;
    for order = [1 2 3]
        y = broaden_binned_centers(x, 0, 1, order, gamma, ...
            'legacy_lorentzian');
        m = profile_fwhm(x, y);
        expected = 2*gamma/order;
        assert(abs(m.fwhm - expected) < 2e-3, ...
            'Legacy Lorentzian width no longer scales as 2*gamma/order.');
    end
    fprintf('  PASS: legacy broadening widths\n');
end
