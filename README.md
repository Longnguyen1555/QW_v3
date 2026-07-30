# QW_MOAP_DeltaDoped

MATLAB project for nonlinear multi-photon optical absorption power in a
Si delta-doped anharmonic GaAs quantum well with:

- self-consistent 1D Schrodinger-Poisson solution;
- perpendicular electric field and in-plane magnetic field;
- electron-LO-phonon scattering with Debye screening;
- electron-piezoelectric-phonon scattering with Debye screening;
- 1PA, 2PA and 3PA;
- phonon emission and absorption separated;
- spectra versus photon energy;
- parameter sweeps versus B, T, Lz, E, doping, alpha or U0;
- FWHM/HWHM extraction without the Signal Processing Toolbox.

## Run

1. Open this folder in MATLAB.
2. Run `Main_quick_test.m` first.
3. Run `Main.m`.
4. Edit the USER SETTINGS block in `Main.m`.

The results are saved in `results/`.

## Main physical equations encoded

1. Anharmonic potential

   `U_A(z) = U0 (z/Lz)^2 [alpha (z/Lz)^6 - 1]`

2. Self-consistent Poisson equation

   `d^2 U_H/dz^2 = e^2/(eps_s eps0) [D(z)-n(z)]`

3. Effective 1D Hamiltonian

   `H = -hbar^2/(2m*) d^2/dz^2
        + (hbar kx + e B z)^2/(2m*)
        + U_A + U_H - e E z`

4. Dipole matrix

   `A_mn = integral psi_m^*(z) z psi_n(z) dz`

5. Electron-phonon form factor

   `I_mn(qz) = integral psi_m^*(z) exp(i qz z) psi_n(z) dz`

6. Optical coupling (screened)

   `V |C_q^op|^2 =
      e^2 hbar omega_LO/(2 eps0)
      (1/eps_inf - 1/eps_s) q^2/(q^2+qd^2)^2`

7. Piezoelectric coupling (screened)

   `V |C_q^pi|^2 =
      kappa^2 hbar e^2 s/(2 eps_s eps0)
      q^3/(q^2+qd^2)^2`

8. Lorentzian replacement of the Dirac delta and numerical cylindrical
   q integration.

## Important modeling choice

The source draft contains copied Landau-level notation even though its
Hamiltonian uses an in-plane magnetic field and solves the z-dependent
states directly. The code therefore does not add a separate Landau
energy ladder. It uses the computed subband energies and wavefunctions.
The Eq. (12) Landau form factor is retained as an optional in-plane
momentum factor because it is needed for article-[6]-style convergent
multi-photon spectra. At B=0 it is automatically replaced by a thermal
recoil cutoff.

## Absolute scale

Raw spectra are saved. The plot is normalized by the global maximum by
default because normalization area/volume conventions in the supplied
draft are not fully specified. No artificial multiplicative `1e4`,
`1000`, or similar linewidth/intensity factors are used.

## Tests

Run:

```matlab
run_all_tests
```

from the `tests` folder or after adding the project recursively to the
MATLAB path.

## Debye screening default

The default screening density is the fixed value `1e18 cm^-3` listed in
the supplied calculation note. Switch to `screening_density_mode =
'from_sheet'` only when a sheet-to-volume conversion is intended.
