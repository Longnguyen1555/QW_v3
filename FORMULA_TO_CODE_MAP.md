# Formula-to-code map

| Physical quantity / equation | MATLAB implementation |
|---|---|
| Anharmonic potential `U_A(z)=U0(z/Lz)^2[alpha(z/Lz)^6-1]` | `src/sp/anharmonic_potential.m` |
| Finite-difference Hamiltonian | `src/sp/build_sp_hamiltonian.m` |
| Charge-neutrality Fermi level | `src/sp/solve_fermi_level.m` |
| Electron density `n(z)` | `src/sp/electron_density_from_subbands.m` |
| Poisson equation for `U_H` | `src/sp/solve_poisson_dirichlet.m` |
| Self-consistent mixing and convergence | `src/sp/solve_schrodinger_poisson.m` |
| Dipole matrix `A_mn` | `src/oap/compute_transition_data.m` |
| Form factor `I_mn(qz)` | `src/oap/compute_transition_data.m` |
| Debye wave number `q_d` | `src/oap/debye_wavevector.m` |
| Optical screened coupling in Eq. (2) | `src/oap/phonon_coupling_density.m` |
| Piezoelectric screened coupling in Eq. (18) | `src/oap/phonon_coupling_density.m` |
| Optical/piezo Bose-weighted terms | `src/oap/direct_phonon_terms.m` |
| Exact positive-k Dirac root and Fermi factor | `src/oap/direct_k_delta_weight.m` |
| Direct Eq. (2)/(18) q-space integral | `src/oap/compute_moap_direct.m` |
| Lorentzian broadening (not used by `direct_q_integral`) | `src/oap/broaden_binned_centers.m` |
| Analytical Taylor-Bessel Eq. (5)/(7) | `src/oap/compute_moap_analytical_series.m` |
| Peak position, FWHM and HWHM | `src/oap/profile_fwhm.m` |
| Spectra versus B/T/L/E/Nd/alpha/U0 | `src/core/run_parameter_sweep.m` |
| Article-[6]-style plots | `src/plot/*.m` |
