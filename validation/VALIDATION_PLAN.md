# Validation plan

The project includes four MATLAB test groups:

1. wavefunction normalization;
2. charge neutrality;
3. direct Eq. (2)/(18) physics: exact k root/domain, Fermi dependence,
   recoil, Bose conventions, acoustic piezo dispersion, photon-order
   validation, gamma independence, q approximation, finite output, and
   two-grid convergence;
4. finite, positive optical and piezoelectric end-to-end spectra.

Additional research-grade checks:

- compare quick, standard and high_accuracy q grids and cutoffs;
- double q cutoffs while holding grid spacing approximately fixed;
- verify the dipole selection rule in the symmetric zero-field limit;
- compare direct_q_integral with analytical_series only after confirming
  convergence in `s_max` and `eta_max`;
- reproduce the qualitative reference-[6] trends:
  1PA peak to the right of 2PA,
  optical emission peak to the right of optical absorption,
  temperature changes intensity more strongly than optical peak position.
