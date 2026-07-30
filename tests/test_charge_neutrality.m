function test_charge_neutrality()
    cfg = apply_numerical_profile(default_config(),'quick');
    sp = solve_schrodinger_poisson(cfg);
    rel = abs(sp.charge_integral_m2-sp.Nd_sheet_m2) / ...
          max(sp.Nd_sheet_m2,1);
    assert(rel < 1e-6, 'Charge neutrality test failed.');
    fprintf('  PASS: charge neutrality\n');
end
