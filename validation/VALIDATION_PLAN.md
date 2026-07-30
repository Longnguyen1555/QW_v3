# Validation plan

The project includes four MATLAB tests:

1. wavefunction normalization;
2. charge neutrality;
3. optical resonance positions
   `(DeltaE +/- hbar*omega_LO)/ell`;
4. finite, positive end-to-end spectrum.

Additional research-grade checks:

- compare quick, standard and high_accuracy profiles;
- double q cutoffs while holding grid spacing approximately fixed;
- verify the dipole selection rule in the symmetric zero-field limit;
- compare direct_q_integral with analytical_series only after confirming
  convergence in `s_max` and `eta_max`;
- reproduce the qualitative reference-[6] trends:
  1PA peak to the right of 2PA and 3PA,
  optical emission peak to the right of optical absorption,
  temperature changes intensity more strongly than optical peak position.
