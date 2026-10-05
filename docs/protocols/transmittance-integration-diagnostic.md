# Protocol: is the remaining longwave gap the way transmittance is integrated?

**Frozen before the calculation.** Results are appended below a line at the
end; nothing above it changes afterwards. This is a diagnostic. It changes no
module.

## Why

`water-vapour-longwave.md` attempt 2 fitted Chou's band 1, 2 and 8 water tables
to 0.0013-0.0017 in transmittance. Even so, it moved the MLS outgoing longwave
further from line-by-line, from 287.80 to 286.48 against 293.99.

**The hypothesis left there, unmeasured, is the integration.**

- **Chou's way.** He integrates band-mean transmittances between pairs of
  levels, each path at its Curtis-Godson effective pressure and temperature
  (eq. 4.13-4.14).
- **This model's way.** The solver treats each bin layer by layer, at each
  layer's own pressure. That assumes the k-terms keep their order from layer
  to layer.

## The calculation

The three bands where water is the only absorber and Chou's tables exist:

- band 1, 0-340 cm^-1;
- band 2, 340-540 cm^-1;
- band 8, 1380-1900 cm^-1.

The atmosphere is MLS as in V1: 1 km layers to 1 hPa, and Chou's 4e-6 g/g of
water above the tropopause. For each band the TOA upward flux is computed
three ways.

- **A, this model's solver.** Per bin, layer by layer, with the fitted
  pressure-tabulated k-distribution of attempt 2.
- **B, Chou's integration, the same coefficients.**
  - **The flux:** F_top = B(T_s) T(top, sfc) + sum over layers of B(T_l)
    [T(top, upper edge of l) - T(top, lower edge of l)].
  - **The transmittance:** T = sum_n g_n exp(-1.66 k_n(p_eff) w h(theta_eff))
    over the path from each interface to the top.
  - **The path:** amount w, with p_eff and theta_eff its amount-weighted means.
  - **The Planck flux:** B is the band's Planck flux at each temperature.
- **C, Chou's integration with his tables.** As B, but the transmittance is
  taken from the tables directly, a + b dtheta + c dtheta^2, interpolated
  bilinearly in log10 w and log10 p_eff.

**Reported for every band:** A, B, C, line-by-line, and the same at 4x layer
refinement.

## The verdict

- **X = B - A,** summed over the three bands, at 1 km: the integration method
  alone, with identical coefficients.
- **Y = C - line-by-line,** summed: whether Chou's method, as implemented here,
  reproduces his own result.

| finding | means |
|---|---|
| \|Y\| <= 1.5 and \|X\| >= 1.5, X closing the gap | the integration carries the gap. The fix is coefficients that hold their order across layers (k-distributions derived from line-by-line spectra, such as RRTMG's), not a different fit to Chou's tables |
| \|Y\| <= 1.5 and \|X\| < 1.5 | not the integration. The gap is elsewhere, and this is reported |
| \|Y\| > 1.5 | the reference implementation of Chou's method is itself off, and nothing is concluded about A |

The sum of the three bands' line-by-line values is 34.56 + 61.01 + 7.56 =
103.13 W m^-2.

---

## Results (appended; nothing above this line changed)

**A unit error, found by reading the first output.** The first run passed B's
effective pressure, already in Pa, multiplied by 100 again. It gave B = 78.77
against C = 100.22, which with identical coefficients could only be a defect.
The first run's X was -19.78 and Y was -2.91; Y does not involve B and is
unchanged. The corrected run is the one reported.

MLS, W m^-2:

| band | A: solver | B: Chou's integration, same k | C: Chou's integration, his tables | line-by-line |
|---|---:|---:|---:|---:|
| 1 | 33.44 | 33.79 | 33.79 | 34.56 |
| 2 | 58.33 | 59.39 | 59.40 | 61.01 |
| 8 | 6.78 | 7.01 | 7.03 | 7.56 |
| **sum, 1 km layers** | **98.55** | **100.19** | **100.22** | **103.13** |
| sum, 4x refined | 98.65 | 100.30 | 100.33 | |

- **X = B - A = +1.64.** The solver's layer-by-layer treatment loses 1.6
  W m^-2 across these bands relative to Chou's path integration, with
  identical coefficients. Refinement does not change it.
- **Y = C - line-by-line = -2.91** (-2.80 refined).

**Verdict as frozen: |Y| > 1.5.** The reference implementation of Chou's method
is itself 2.9 W m^-2 short of his line-by-line. It is also 2.3 short of his own
HIGH-option result, 34.40 + 60.54 + 7.61 = 102.55. So nothing is concluded
about A.

**Reported.**

- **B and C agree to 0.03.** The fitted coefficients reproduce Chou's tables
  along real, inhomogeneous paths, not only on the homogeneous ones they were
  fitted to.
- **The benchmark is not quite like for like.** Chou's section 10 says his
  line-by-line used the McClatchey et al. (1972) mid-latitude summer, in 75
  layers (25 hPa below 100 hPa). This uses AFGL 1986, Table 1b. His water
  above the tropopause, 4e-6 g/g, matches AFGL's 3-8 ppmv between 179 and
  100 hPa, so that is not the difference.
- **A shortfall shared by all three.** A, B and C share the atmosphere, and all
  three fall short of line-by-line together. So a common part of the gap is
  plausibly the atmosphere or the remaining details of Chou's integration,
  such as his separate treatment of adjacent layers. That is not measured.
- **How much is attributable,** for these three bands of the -6.2 W m^-2 total
  (attempt 1 MLS):
  - 1.6 to the solver's integration;
  - up to about 2.3 to things common to any implementation run on this
    atmosphere.
