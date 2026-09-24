# Validation plan

The project includes five MATLAB tests:

1. wavefunction normalization;
2. charge neutrality;
3. optical resonance positions
   `(DeltaE +/- hbar*omega_LO)/ell`;
4. refined Eq. (32) versus Eq. (28), independent adaptive quadrature,
   specific invalid-series errors and transition scaling;
5. finite, positive end-to-end spectrum.

Additional research-grade checks:

- compare quick, standard and high_accuracy profiles;
- double q cutoffs while holding grid spacing approximately fixed;
- verify the dipole selection rule in the symmetric zero-field limit;
- refine `s_max` and `v_max` independently for Eq. (32), and `Nqperp`
  for its Eq. (28) reference; numerical agreement does not resolve the
  formula differences documented in `ANALYTICAL_SERIES_REFINED_AUDIT.md`;
- reproduce the qualitative reference-[6] trends:
  1PA peak to the right of 2PA and 3PA,
  optical emission peak to the right of optical absorption,
  temperature changes intensity more strongly than optical peak position.
