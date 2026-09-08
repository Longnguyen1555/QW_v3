# Direct q-integral physics audit

## Scope

`cfg.oap.model = 'direct_q_integral'` dispatches through
`compute_moap_spectrum.m` to `compute_moap_direct.m`.  It is distinct
from `direct_t`, which dispatches to `compute_moap_pdf_direct.m`.
The Schrödinger-Poisson and B-field Hamiltonian are outside this patch.

## Eq. (2): LO optical phonons

The direct path uses constant `hw0`, the screened coupling

`V|Cq|^2 = e^2*hw0/(2*eps0)*(1/kappa_inf-1/kappa_0)
            *q_perp^2/(q_perp^2+qd^2)^2`,

and Bose factors `N0+1` for emission and `N0` for absorption.  For each
photon energy and `ell=1,2`, the laser factor is

`(a0*q_perp)^(2*ell)/(2^(2*ell)*factorial(ell)^2)`,

with `a0=eE0/(m*Omega^2)`.

## Eq. (18): piezoelectric acoustic phonons

The screened coupling changes its numerator to `q_perp^3`.  The written
source equation reuses `N0` and `omega0` from the LO section despite
introducing sound speed `c_s`.  The physical direct implementation makes
that correction explicitly:

`hw_q=hbar*c_s*q_perp`, `N_q=1/expm1(hw_q/kBT)`.

At `q_perp=0`, `N_q` diverges, but the complete screened product
`C2V*N_q` tends to zero for finite Debye screening; that product limit,
not a fictitious zero Bose population, is used.

## Exact k integral and source approximations

With the source's collinear convention `k dot q = k*q`,

`B(q)=hbar^2*q^2/(2m*)+DeltaE+zeta*hw_q-ell*Eph`,

`kstar=-m*B/(hbar^2*q)`.

The radial domain is `k in [0,Inf)`, so a contribution exists only for
`kstar>=0` (equivalently `B<=0` for `q>0`).  The exact delta integral is

`m*kstar/(hbar^2*q)*f(E_initial+hbar^2*kstar^2/(2m*)-EF)`.

The source omits the required positive-root domain statement.  The code
enforces it.  It deliberately does not add a general angular integral.

The approximation `q approximately q_perp` is used in coupling, acoustic
dispersion, and laser dressing.  `q_z` remains only in
`|I_n'n(q_z)|^2`.  Because that grid covers `[0,Inf)`, a factor of two
restores the negative `q_z` half-space.

## Line shape and normalization

The direct path does not call `broaden_binned_centers`; configured gamma
values do not set its line shape.  FWHM is extracted after the numerical
q integration.

The normalization volume in `|Cq|^2` cancels the 3-D q-state density.  The
in-plane electron sum contributes `S/(2*pi)` before the radial k integral,
where `S=cfg.oap.normalization_area_m2`.  The repository retains the
source electromagnetic/golden-rule prefactor convention.  That convention
is not dimensionally closed as an all-SI absolute-W expression in the
provided derivation, so raw spectra are explicitly labeled source-normalized
and no absolute-W claim is made.

In SI bookkeeping, `V|Cq|^2` has units `J^2 m^3`; the q-state integral
removes `m^3`; and `S/(2*pi)` times the exact k/delta integral has units
`1/J`.  Before the source electromagnetic factor, the golden-rule factor
therefore produces a rate.  Substituting the configured `E0` in `V/m`
leaves uncancelled field units, which is the specific unresolved reason that
the result is not labeled as absolute watts.
