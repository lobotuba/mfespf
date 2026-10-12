# Protocol: the longwave source function, correlated-k emission and linear-in-tau

**Frozen before any code, check or run.** Results are appended below a line
at the end; nothing above it changes afterwards.

## Why

`rrtmg-longwave.md` P1 failed on three of four atmospheres. Its
decomposition showed that the gas optics are RRTMG's exactly (P0) and that
the gap lies in how this model's solver makes the layers emit. Two of its
three parts are physics this solver does not have:

- **(p), correlated-k emission.** A g-point is not a wavelength interval. It
  is the set of wavenumbers throughout its band whose absorption falls in one
  range. Its emission is therefore its Planck fraction times the band's
  Planck flux (Lacis & Oinas 1991, *JGR* 96, 9027; Mlawer et al. 1997, *JGR*
  102, 16663).
  - This model places each g-point on a wavelength sub-interval and lets it
    emit that sub-interval's Planck flux. Its share of the band then shifts
    with temperature, which a g-point's does not.
  - P1 measured -0.9 to -2.4 W m^-2 of OLR from this.
- **(s), the linear-in-tau source.** A layer's emission is not at one
  temperature. An optically thick layer emits from near its edges, at their
  temperatures, and a thin one from its whole depth.
  - RRTMG lets the Planck function vary linearly in optical depth between the
    level temperatures (Clough, Iacono & Moncet 1992, *JGR* 97, 15761; RRTMG
    `rtrnmr`). This model's layers emit at their mean temperature.
  - P1 measured -7.1 W m^-2 of the MLS surface downward flux from this.

The third part, (d), is RRTMG's band-dependent diffusivity: AER's fitted
`secdiff`, 1.5-1.8 by band and column water. It is not built here, and its
measured size is the residual this protocol predicts.

## The physics

**(p) Correlated-k emission, `k_distribution_source`.** It applies only to
bins that hold an RRTMG g-point.

- **The layers.** In layer i, the bin holding g-point g of band b emits
  `F[i, g] x B_b(T)`.
  - F[i, g] is RRTMG's own Planck fraction for that layer, from `taumol`.
    It is the same array whose fixed version sized the grid.
  - B_b is the band's Planck flux over the band as it stands on the grid,
    the sum of its bins' Planck fluxes. Band 15 is cut at 2500 cm^-1.
- **The surface** emits `F[0, g] x B_b(T_s)` per bin, times its emissivity.
  That is `rtrnmr`'s choice: the lowest layer's fractions.
- **Every other bin** keeps its own Planck flux.
- **The fractions must close.** Each band's surface fractions must sum to one
  to 1e-6, or the solver refuses. The surface module emits its own greybody
  flux, so the radiation's surface emission must carry the same total; it
  does exactly when the fractions close.
  - RRTMG's lower-atmosphere fractions close on all four of AER's
    atmospheres, checked before this protocol was frozen.
  - Above laytrop, the fractions of bands 12 and 15 are zero, and so are
    those bands' optical depths, so no energy is created.

**(s) The linear-in-tau source, `linear_in_tau_source`.** It applies to every
bin.

- **The source.** A layer of emissivity e emits upward
  `S_up = e [B_lay + f(x)(B_top - B_lay)]`, and downward the same with
  B_bottom.
  - x = D tau_a, where tau_a is the layer's absorption optical depth: gas
    plus cloud absorption, not scattering.
  - D is the diffusivity factor, unchanged.
  - B_top and B_bottom are the Planck fluxes at the interface temperatures.
    Under (p) they are weighted by the layer's own fractions, as `rtrnmr`
    does.
- **The weighting function**, from RRTMG `rrtmg_lw_init`:
  - f(x) = 1 - 2 (1/x - e^-x / (1 - e^-x)) for x > 0.06;
  - f(x) = x/6 below.

  It is evaluated directly rather than tabulated.
- **Its limits:**
  - a transparent layer, f → 0, emits at its mean;
  - an opaque layer, f → 1, emits at its edge temperature;
  - an isothermal layer, B_top = B_lay = B_bottom, emits exactly what the
    current solver gives.
- **Where it is approximate.** In a layer that does not scatter, this is
  `rtrnmr`'s source exactly. In a scattering layer, such as a cloud, the
  weighting is applied to the two-stream emissivity e = 1 - r - t. That is an
  approximation and is declared as one. RRTMG's own longwave does not scatter.
- **The adding method** carries the two sources separately: S_up leaves the
  layer's top and S_down its bottom. The recurrences are the present ones,
  with S_down in the reflected term and S_up in the transmitted.

**The interface temperatures.** These follow CAM's RRTMG driver (ESCOMP/CAM
`src/physics/rrtmg/rrtmg_state.F90`, "following method used in radtpl"):

- **The top interface:** the top layer's temperature.
- **The interior interfaces:** linear in ln p between the adjacent layer
  temperatures.
- **The bottom interface:** CAM's (LW_up/σ)^1/4. For this model that is the
  surface temperature, whose longwave emissivity is one.

**Both switches default off.** With both off, the solver is unchanged to the
bit.

## The tests

**Unit tests, all must pass before P2:**

- **T1, identity:** both off, every flux is identical to the present solver's.
- **T2, isothermal:** with the interface temperatures equal to the layer
  temperatures, the linear source equals the isothermal one to 1e-12
  relative.
- **T3, the limits of f:** f(x)/x → 1/6 as x → 0, f → 1 as x → ∞, and f is
  continuous at 0.06 to 1e-4.
- **T4, an opaque layer:** a single layer with D tau = 50, a 280 K top and a
  300 K bottom emits upward B(280) and downward B(300), to 1e-6 relative.
- **T5, fixed fractions close:** with layer-independent fractions, an
  isothermal column over a black surface at the same temperature gives each
  band's OLR equal to B_b(T), to 1e-9 relative.
- **T6, refusals:**
  - `k_distribution_source` without `rrtmg_longwave`;
  - fractions that do not close;
  - interface temperatures of the wrong length.
- **T7, the module:** with both switches on, the column's energy residual
  stays at round-off.

**P2, decisive: does this solver now reproduce RRTMG's transfer?** It is P1
exactly:

- AER's four atmospheres, the grid to 2600 cm^-1;
- the reference is AER's total minus band 16;
- the thresholds are unchanged: OLR within 1.5 W m^-2 and surface downward
  flux within 2.0.

Both switches are on. The interface temperatures are AER's own level
temperatures `tz`, and the surface is AER's `tbound`.

- **Passes if** all four atmospheres meet both limits.
- **Predicted before the code,** from P1's decomposition: the residual is
  (d) alone, W m^-2:
  - MLS -0.60 / +0.24;
  - SAW -0.16 / +0.49;
  - TROP -0.97 / +1.12;
  - MLW -0.23 / +0.07.

  The tolerance on that prediction is 0.15 each, for the direct functions in
  place of `rtrnmr`'s tables. P2 should pass on all four.
- **Reported:**
  - each switch alone;
  - the CAM interface convention in place of AER's `tz`, which measures the
    convention's cost;
  - RRTMG's fixed fractions in place of the per-layer ones.

**V1, V2, V5, re-judged.** These are `rrtmg-longwave.md`'s, exactly, with
both switches on. The interface temperatures are the AFGL level values. The
limits are unchanged, and so is the caveat recorded there: they measure the
disagreement between RRTMG and Chou's line-by-line as much as this model.

**R4, the global column.** R3's configuration exactly, plus both switches.

- **The run:** 200000 steps, with the state saved.
- **G0-G4** as frozen in `global-column-chou-longwave.md`. G4 is against
  R0's +1.59 K.
- **Predicted before the run:** G4 still fails. This is not a solver
  problem: the attribution put 39.7 W m^-2 on the excess water, and nothing
  here dries the column. The surface stays within 3 K of R3's 302.65 K.
  - (p) lets out about 2 W m^-2 more at the top for a given state.
  - (s) sends more to the surface from thick, moist lower layers.
  - The sign of the net is not predicted.
- **Reported:** everything R3 reported, against R3.

## Order and outcomes

**The order:** the unit tests, then P2, then V1, V2 and V5, then R4.

- **If P2 fails,** the cause is found and recorded. A defect in transcribing
  a stated formula is corrected and recorded. Nothing here may be tuned:
  - not the fractions;
  - not f;
  - not the interface convention;
  - not D.

  R4 runs only on a P2 pass.
- **R4's result decides nothing by itself.** Adoption of either switch is the
  user's decision.
- **Then:** the free-tropospheric drying mechanism, under its own protocol,
  starts from R4's configuration whatever R4's result.

---

## Erratum, before any result (nothing above this line changed)

**The closure tolerance was set from a check printed to four decimals.** The
protocol says each band's surface fractions must sum to one within 1e-6.

**Measured on AER's four atmospheres:** RRTMG's own lower-atmosphere
fractions miss one by 1.5e-6 to 1.72e-5, depending on the band. That is the
rounding of the published tables, the same on every atmosphere. As frozen,
the solver refuses RRTMG's own fractions, and (p) cannot run at all.

**The two readings:**

1. **Loosen the tolerance to admit the rounding.** This weakens a check, so
   it is not taken.
2. **Keep the check at 1e-6 and divide the surface fractions by their band
   sum,** taken.
   - This restores the defining property of a g-point's fractions, which sum
     to one over the band, where the tables lose it to rounding.
   - `rrtmg_fixed_fractions` already does the same per band
     (`tools/gen_rrtmg_fractions.py`).
   - It changes the surface's emission by at most 1.72e-5 of each band's
     flux.
   - The layer fractions are used as published. Heating is flux divergence,
     so the column's energy closes identically whatever they sum to.

**What is recorded with P2:** the size in W m^-2 of the normalisation, as
the difference between the normalised and published surface emission.

## Second erratum, T4, before any result is recorded (nothing above this line changed)

**T4 contradicts the formula frozen with it.**

- **What T4 says:** a layer at D tau = 50 emits upward B(280) and downward
  B(300) to 1e-6.
- **What the frozen formula gives:**
  - f(x) = 1 - 2 (1/x - e^-x / (1 - e^-x)), so f(50) = 0.96, not 1;
  - the upward emission is B_lay + 0.96 (B_top - B_lay);
  - the emission approaches the edge's only as 1/x.

  That is the physics: an opaque layer emits from about one optical depth
  inside its edge, so the departure is (B_top - B_bottom)/x.
- **The first run:** the code gives the formula's value. T4 as written fails
  by 0.6 percent at the spectral peak.

**The two readings:**

1. **Keep T4 as frozen.** It fails, and it judges the code against a value
   the frozen physics does not predict.
2. **Test the frozen formula where it says it holds,** taken:
   - at D tau = 50, upward and downward equal
     B_lay + f(50)(B_edge - B_lay) per bin, to 1e-6 relative;
   - at D tau = 1e9, the limit, the total upward and downward equal the
     edge's Planck flux, B(280) and B(300), to 1e-6 relative.

   Both are exact consequences of the formula. Neither loosens a tolerance.

## Unit tests, P2 and V result (appended; nothing above this line changed)

**The unit tests pass,** under both errata:

- T1-T7;
- CAM's interface convention against a hand calculation;
- P2 as a permanent test, `test_the_solver_reproduces_rrtmg_s_transfer`.

**The equation records** are `EQ-CORRELATED-K-EMISSION` and
`EQ-LINEAR-IN-TAU-SOURCE`. The switches are `k_distribution_source` and
`linear_in_tau_source` on the radiation module and `reference_column`. Both
are off by default.

**P2 passes on all four atmospheres, and the residual is the predicted one.**
All values in W m^-2, OLR / surface down, against AER's bands 1-15:

| atmosphere | P1, neither | **P2, both** | predicted, (d) alone | P2 |
|---|---:|---:|---:|---|
| MLS | -2.27 / -6.00 | **-0.60 / +0.24** | -0.60 / +0.24 | pass |
| SAW | -1.06 / +0.75 | **-0.16 / +0.49** | -0.16 / +0.49 | pass |
| TROP | -3.27 / +0.47 | **-0.97 / +1.12** | -0.97 / +1.12 | pass |
| MLW | -1.62 / +0.13 | **-0.23 / +0.07** | -0.23 / +0.07 | pass |

Every residual is within 0.01 W m^-2 of the prediction made from P1's
decomposition before the code. What remains is RRTMG's band-dependent
diffusivity, which is not built.

**Reported as frozen:**

| atmosphere | (p) only | (s) only | both, CAM's interfaces in place of AER's `tz` | both, fixed fractions in place of per-layer |
|---|---:|---:|---:|---:|
| MLS | -0.07 / -6.91 | -2.80 / +1.29 | -0.61 / +0.24 | +0.04 / +0.20 |
| SAW | -0.12 / +0.63 | -1.10 / +0.61 | -0.16 / +0.49 | +0.10 / +0.45 |
| TROP | -0.85 / -0.48 | -3.38 / +2.11 | -0.97 / +1.12 | -0.19 / +1.05 |
| MLW | -0.18 / -0.39 | -1.68 / +0.60 | -0.23 / +0.07 | +0.11 / +0.03 |

- **(p) carries the OLR and (s) the MLS surface,** as P1 found. Each alone
  leaves the other's error in place.
- **CAM's interface convention costs at most 0.01.** AER's surface level
  temperature equals its `tbound` on all four atmospheres, so that is
  expected.
- **The erratum's normalisation** changes the surface emission by -0.0004 to
  -0.0009 W m^-2.

**V1, V2 and V5, re-judged with both options.** The interface temperatures
are the AFGL levels. All values in W m^-2:

| atmosphere | V1: OLR - Chou line-by-line | V2: surface down - Chou line-by-line | V5 |
|---|---:|---:|---|
| MLS | 283.33 - 293.99 = **-10.66**, fail | 348.14 - 339.59 = **+8.55**, fail | -12.43 C/day at 1.1 hPa, pass |
| SAW | 200.60 - 204.37 = **-3.77**, fail | 171.01 - 162.01 = **+9.00**, fail | (-27.7) |

**MLS's V2 passed before (+2.93) and fails now.** The earlier pass was a
cancellation:

- the solver's own -6.0 on the surface flux, which P2 removes;
- against RRTMG's +8.5 disagreement with Chou's line-by-line, which
  `rrtmg-longwave.md` measured through AER's own transfer (348.09 on this
  column).

The model's solver now gives 348.14 there. What the V tests measure now is
the disagreement between the two references, as recorded.

**As frozen, R4 follows (P2 passed).**
