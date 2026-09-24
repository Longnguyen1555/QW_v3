# Refined analytical-series formula map

Source: `Refined Tính_toán_chi_tiết.pdf`.

| PDF expression | Implementation |
|---|---|
| Eq. (20), optical integral in its second printed form | `compute_refined_T_pdf`, `parent_integral` |
| Eq. (22), LO power | `compute_moap_analytical_series`, `P0_optical` and LO Bose factors |
| Eq. (24), piezoelectric `D` | `Dpiezo = DeltaE - ell*Ephot` |
| Eq. (26)-(30), piezoelectric power | `P0_piezo` and the combined acoustic contribution |
| Page 7 Taylor/Fermi sums | `bessel_sum`, `s = 0:s_max`, `v = 1:v_max` |
| Eq. (32), Bessel-K result | `bessel_sum` with scaled `besselk` |

## Internal inconsistencies resolved

The PDF contains three compact-form typographical inconsistencies. The code
uses the algebraically consistent result obtained from the preceding steps:

1. Eq. (22) prints `E0^2`, whereas combining Eqs. (2), (4), (5), (7),
   (11), and (18) gives `E0^4`.
2. Eq. (27) prints `Omega^-1`, whereas Eq. (5) and the collected prefactor
   on page 5 give `Omega^-2`.
3. The double sum on page 7 prints `(v+1)`, whereas the immediately preceding
   Taylor expansion of `(1+gamma^2)^-2` gives `(s+1)`.

No Lorentzian broadening, clipping, or arbitrary intensity rescaling is used.

Another source inconsistency matters for negative D: Eqs. (19)/(25) contain
`abs(Eq + D)`, while Eqs. (20)/(28) and Section 4 use `Eq + abs(D)` with
`Eq = hbar^2*q^2/(2*m)`. These are NOT equal when D < 0. Both numerical
implementations deliberately retain the latter expression to match Eq. (32).
Agreement between them does not validate that earlier change of formula.

## Review of the series-domain error and simplified implementation

`cfg.oap.series.evaluation` selects:

- `auto` (default): compute Eq. (28); try Eq. (32) only when `E_initial > EF`
  and `D ~= 0`, and use it only if it is finite, nonnegative and within
  relative error `1e-3` of the Eq. (28) quadrature. Otherwise retain Eq. (28).
- `series`: require the same checks, return the actual Eq. (32) sum, and
  fail with a specific error when a check fails. Never silently substitute
  quadrature for the returned series result. Quadrature is still computed
  as the independent numerical reference for this check.
- `quadrature`: use Eq. (28) only; `s_max` and `v_max` do not affect it.

The old `Eqd/abs(D) <= 0.25` gate was a heuristic, not a source-derived
convergence condition. It has been removed. The old comparison of two sums
with both indices reduced together could hide cancellation; it has also
been removed. The `1e-3` agreement tolerance is a numerical acceptance test,
not a physical cutoff or a proof of convergence of the infinite series.

The Taylor series in `gamma=qd/q` converges only for `abs(gamma) < 1`.
The substitution does not enforce this on `q in (0,Inf)`. In particular,
the termwise integrated screening series must not be assumed globally
convergent merely because a finite truncation is positive or changes little.
The PDF's stated sufficient Fermi condition `E_initial > EF` is enforced
for series evaluation. Increasing truncation orders cannot repair a violation.

For the current saved 90 K case, `E1=-298.751965 meV`, `EF=-119.093455 meV`,
and `exp((EF-E1)/kBT)=1.1493e10`. Forced series therefore correctly stops;
`auto` evaluates Eq. (28) for this case. No electronic parameters are altered.

Quadrature uses `max(96,2*Nqperp)` Gauss-Legendre nodes on the transformed
infinite domain. Refine `Nqperp` for an accuracy study, especially at low T.
The cancellation-resistant square `(Eq+D)^2/(4*Eq)` evaluates Eq. (30).
This is algebraically identical to `Eq/4 + D/2 + D^2/(4*Eq)`.

Unused `eta_max`, `piezo_effective_q_mode` and the Fermi-check bypass flag
were removed from the default configuration; extra fields in old saved
configurations are harmless and ignored. Diagnostics now record just the
counts of returned series and quadrature values.

The acoustic Eq. (26) result remains in the legacy `emission_raw` slot
with `absorption_raw=0`, as recorded in metadata. It is a combined acoustic
contribution, not a separate emission/absorption prediction. Existing
generic plotting/CSV labels still use the legacy slot names.

## Verification of the simplified implementation

- The two kernel files were reduced from 407 to 186 lines in total.
- All five MATLAB tests pass; Code Analyzer reports zero messages for the
  changed kernels, configuration and targeted test.
- The targeted test covers both signs of D, both photon orders, an
  independent dimensionless adaptive integral, the three specific series
  errors, dipole-square scaling and addition of multiple transitions.
- On the same saved electronic input, the quadrature spectrum differs
  from the pre-refactor result by `3.07e-16` relative to its maximum.
- A fresh standard SP calculation at 90 K with `Lz=5 nm`, `Nd=1e13 cm^-2`,
  `B=0`, and `2:0.01:250 meV` completes in auto mode with finite,
  nonnegative values. All 148806 T evaluations use quadrature because the
  initial subband is below EF. The spectrum calculation took 2.09 s.
- Refining the quadrature from 180 to 360 nodes changes the spectral area
  by `2.94e-6` relative and the maximum absolute difference divided by the
  refined maximum is `1.63e-6`. This is a single 90 K validation, not a
  convergence claim for the entire temperature sweep.
