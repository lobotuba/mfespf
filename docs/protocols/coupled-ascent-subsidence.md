# Protocol: a dry subsiding branch, from a fixed contrast to prognostic surfaces (S1, S2)

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards. The design is the
user's (2026-10-17), recorded here with its implementation fixed.

**It follows from:**

- `tropical-runaway-trajectory.md`. N0's own trajectory confirmed the
  radiative-moisture attribution: no free-tropospheric drying, RH about
  0.85-0.95, too little OLR, warming toward a hot equilibrium near 320 K.
- `spectral-wtg-circulation.md`. Two identical columns never make a
  temperature difference, so the relaxed WTG balance never circulates.
- `ocean-heat-transport-circulation.md`. W1-W4 ran away together before the
  ocean's contrast could act, and evaporation erased the surface contrast.

## The earlier decision, resolved (recorded as the user worded it)

> The earlier ocean-transport decision applies to the final interactive
> architecture: the sustainable surface-temperature contrast must be predicted
> from the surface-energy budgets and ocean heat transport. S1 does not alter
> that decision. S1 is a diagnostic boundary-condition experiment that imposes
> the observed temperature contrast to test whether the atmospheric coupling
> contains the predicted subsidence-drying pathway. No S1 equilibrium is
> adopted as a coupled climate solution. S2 restores prognostic surfaces and
> tests whether the pathway survives under the frozen observed ocean
> heat-transport contrast.

**Fixed non-uniform surface temperatures are the established WTG way** to
isolate how convection and the parameterized circulation respond to a
contrast (Daleu et al. 2016, *JAMES*, part 2 of the GASS WTG
intercomparison). That study also found single-column responses less regular
than cloud-resolving ones, which is why S1 is treated as a diagnostic.

**Constraints, both stages:**

- current default column physics: no bulk plume, no CAPE-definition
  experiment;
- no humidity relaxation;
- no imposed subsidence;
- no fitted ocean transport;
- no default changed;
- F0, interactive-radiation F0 and P1-P3 stay blocked.

## The observations used (from `ocean-heat-transport-circulation.md`)

±20°, NCEP/NCAR R1 1991-2020 annual means. The ocean heat comes from OAFlux
and NCEP.

| quantity | rising | sinking | band |
|---|---:|---:|---:|
| area fraction | **0.558** | 0.442 | |
| skin temperature | **300.26 K** | **298.79 K** | 299.61 K |
| omega500, hPa/day | −20.2 | **+14.3** | |
| RH 300-700 | 0.50 | **0.35** | **0.43** |
| precipitable water 1000-300, kg m⁻² | 42.3 | 32.7 | **38.0** |
| T500, rising − sinking | | | **+0.65 K** |
| ocean heat delivered, W m⁻² (OAFlux) | +5.68 | −7.71 | contrast **13.39** |
| ocean heat delivered, W m⁻² (NCEP) | +2.03 | −2.64 | contrast **4.67** |

**The sinking region's RH, by level:**

| level, hPa | 1000 | 925 | 850 | 700 | 600 | 500 | 400 | 300 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| RH | 0.75 | 0.69 | 0.57 | 0.39 | 0.32 | 0.28 | 0.28 | 0.50 |

## Code, introduced by this protocol only

**`prescribed_temperature` on `SlabSurfaceModule`,** a float or none. None
is today's code, bit for bit.

- **When set, the surface temperature is held.** At commit, the sum of every
  contribution to the surface temperature this step is cancelled by one
  labelled contribution, `fixed_surface_boundary_transfer`. Its holding
  power is F_hold = −C_s·Σ(dT_s/dt) in W m⁻², positive into the surface.
- **The transfer goes to a boundary reservoir module,
  `fixed_surface_reservoir`, per region.** It is the deep-ocean exchange's
  pattern: the reservoir records the energy cumulatively, and the
  conservation reconciler closes the system including it.
- **The temperature at commit equals the prescribed value to rounding,**
  and the state budget records the cancellation.
- **Published per region:** `fixed_surface_boundary_transfer` (W m⁻²,
  positive into the surface).
- **Passed through the configuration** as `prescribed_surface_temperature`,
  per region of `two_column`.

**Unit tests:**

- option-off identity;
- the held temperature;
- F_hold equal to minus the other contributions;
- conservation closing including the reservoir.

## S1: the drying pathway under a fixed contrast

**The question.** Given the observed tropical area fractions and surface
contrast, does the spectral WTG coupling generate compensating subsidence
that dries the descending free troposphere, before local moistening and
convection erase the dry state? It does not ask whether the model can
predict the contrast. That is S2's question.

**The configuration:**

- `two_column`, `tropical_latitude = 20`;
- ascent fraction 0.558;
- relaxation distance L = 734 km, the ±20° value;
- **surfaces held at 300.26 K (rising) and 298.79 K (sinking);**
- no ocean heat transport;
- W1's column physics: `oht_run.py`'s physics dictionary, master step 3600 s,
  US Standard Atmosphere initial temperatures.

**The length.** 20000 steps, read over the last fifth.

- **Equilibrated if,** in each column over the last fifth, the atmosphere's
  stored-energy rate (from the reconciler's stores) has a mean within ±0.1
  W m⁻², and descent RH 300-700's 50-day means vary by less than 0.02.
- **Otherwise** continued in 20000-step increments, to at most 60000.

### S1 gates

**G1, conservation.** Mass and water close, and atmospheric energy plus the
declared boundary transfers closes: the reconciler passes with a worst
relative residual of 1e-12 or less. The circulation's between-column transfers
cancel, area-weighted.

**G2, the fixed temperature.** Both surfaces stay at their declared
temperature to 1e-9 K, and F_hold accounts for all the holding energy. The
state budget's net change is zero.

**G3, the WTG exchange.** The area-weighted sum of each region's circulation
tendency is zero to rounding, separately for:

- dry static energy;
- vapour;
- total energy (dry static plus latent).

No condensate is exchanged; that is confirmed.

**G4, no local runaway.** Both columns stay inside every module's validity
domain (the validity gate, including RRTMG's ±30 K) for the whole run. Also
bounded over the last fifth, as 50-day means drifting by less than 1 percent:

- the atmospheric temperatures;
- the precipitable water;
- the radiative and adjustment tendencies.

**G5, orientation.** The descent region's omega500 is positive, and the
warm column rises, in at least 90 percent of the daily means of the last
fifth.

**G6, the drying budget.** In the descent column's free troposphere (layers
centred 300-700 hPa), the last fifth's mean of each vapour label, kg m⁻²
day⁻¹, summed over those layers:

- the circulation's vapour tendency;
- `convective_vapour_transport`;
- `vapour_eddy_diffusion`;
- `evaporation_minus_condensation`;
- `brewer_dobson_upwelling`.

**The decisive diagnostic:** whether the circulation's tendency is negative
and survives the local moistening.

**G7, the frozen predictions** (none adjusted; L, σ, the contrast and the
relaxation are never changed to fit them):

| prediction | passes if |
|---|---|
| descent omega500 | 7-29 hPa/day |
| descent RH 300-700 | 0.25-0.45 |
| band RH 300-700, area-weighted | 0.33-0.53 |
| the descent RH profile is C-shaped | its minimum lies between 350 and 650 hPa, below both its 850 and 300 hPa values |
| band precipitable water 1000-300 | 28.5-47.5 kg m⁻² (±25 percent) |
| T500, rising − sinking | +0.33 to +1.30 K (half to double +0.65) |

RH is evaluated at NCEP's levels as `oht_run.py` does.

**Reported, not gated:**

- **the implied surface-transfer contrast,** ΔF = F_hold,rising −
  F_hold,sinking (W m⁻² of each region, positive into the surface). It is
  compared with OAFlux 13.39 and NCEP 4.67, the ocean delivering heat into the
  rising region. It is never used to scale anything;
- **the band mean of F_hold,** σ·F_hold,r + (1 − σ)·F_hold,s, reported
  separately. A large non-zero value means S1 needs net energy added to or
  taken from the whole band, not only redistributed;
- **each column's TOA budget and the applied export;**
- time series of F_hold, omega500 and descent RH.

### S1 verdicts, frozen

| verdict | conditions |
|---|---|
| **Drying pathway supported** | (1) G5 holds; (2) descent omega500 in 7-29 hPa/day; (3) descent RH 300-700 in 0.25-0.45 over the last fifth, with at least 90 percent of daily means ≤ 0.55; (4) band RH 300-700 in 0.33-0.53; (5) in G6 the circulation's tendency is negative and is the largest drying term, larger in magnitude than the negative part of `evaporation_minus_condensation`; (6) G1, G2 and G3 pass; (7) G4 passes. |
| **Local moistening blocker** | Descending motion is sustained (G5). The circulation's free-tropospheric tendency is a measurable drying (≤ −0.1 kg m⁻² day⁻¹). But descent RH 300-700 stays ≥ 0.75, the ceiling's threshold, because convective vapour mixing or another local operator offsets it. Recorded as: *the coupled circulation supplies subsidence, but the default single-column physics prevents the descending branch from retaining its dry state.* |
| **Surface contrast insufficient** | G5 fails, or descent omega500 is below 2 hPa/day over the last fifth. The contrast is not increased within this protocol; that would be a sensitivity study. |
| **Diagnostically successful but energetically unrealistic** | Drying is supported, but ΔF has the opposite sign to the observed delivery (that is, ΔF ≤ 0), or a magnitude below 2.3 W m⁻² (half the NCEP value) or above 26.8 (twice the OAFlux value). |
| **Anything else** | reported as found |

**Only "drying pathway supported", or "diagnostically successful but
energetically unrealistic", opens S2.**

## S2: does the dry circulation maintain itself on prognostic surfaces?

**Run only if S1 opens it.**

**The initial state is S1's final state, directly:**

- each column's temperatures and vapour;
- its surface temperature;
- its water store.

The circulation is diagnosed from the state each step and carries no memory.
Nothing is constructed, nudged or post-processed.

**The transition:**

1. the fixed temperatures are removed;
2. each slab's store is C_s·T_s at S1's final temperature, with C_s = 1e6
   J m⁻² K⁻¹ as in W1;
3. the observed ocean heat-transport contrast is applied: **13.39 W m⁻²
   (OAFlux),** W1's verdict value, through `OceanHeatTransport`;
4. the atmosphere is kept as it is, with nothing reset.

**The energy jump across the transition** is recorded: the total stored
energy at the end of S1 against the start of S2. It must be zero to 1e-12
relative.

**The run.** 20000 steps, read over the last fifth, extended in 20000-step
increments to at most 60000 if the band's net TOA flux is not within ±0.1
W m⁻² of its applied export.

### S2 outcomes, frozen

| outcome | conditions |
|---|---|
| **Self-maintaining coupled state** | Over the last fifth: (1) the rising-minus-sinking surface contrast is positive, between 0.5 and 2 times S1's 1.47 K, and drifts by less than 0.1 K per 50 days; (2) the orientation holds (G5); (3) descent RH 300-700 ≤ 0.55; (4) neither column runs away: the band's mean surface temperature is within 2 K of S1's band value, the band is equilibrated as above, and the validity gate passes; (5) conservation passes. |
| **Collapse to moist runaway** | descent RH 300-700 ≥ 0.75, or descent omega500 < 2 hPa/day, or the surface contrast < 0.5 K, and the band's surface temperature more than 2 K above S1's and rising |
| **Dependence on fixed-boundary support** | The collapse criteria are first met within S2's first 2000 steps. That separates "a dry circulation is possible under prescribed surfaces" from "the interactive system can sustain it". |
| **Anything else** | reported as found |

**Reported:**

- the time series of the surface contrast, omega500, descent RH and the
  band's surface temperature;
- each column's surface energy terms and its ocean heat.

## What this cannot separate

- **The regional partition.** Two regions stand for the ±20° band.
  Variability within each region, and land, are not represented.
- **The clouds are the band's.** The ISCCP tropical cloud field is
  prescribed in both columns. Real clouds differ by regime.
- **S1's verdict is about this coupling.** It is not about the WTG
  literature's other forms (strict, or damped gravity wave).

---

## Erratum to the code section, before any code (nothing above this line changed)

**Reading 1 (above).** The holding energy goes to a boundary reservoir
module, `fixed_surface_reservoir`, per region, in the deep-ocean exchange's
pattern.

**Reading 2, from the coupler's lifecycle, read before writing it.**

- **Why reading 1 cannot work as written.** Every module's tendencies are
  computed and published first, then each owner commits, then transfers and
  budgets are reported. The holding power F_hold is the sum of every
  contribution to the surface temperature, so it is known only at the
  surface's commit. A separate reservoir module would have to read it in the
  tendency pass, before it exists.
- **The established pattern for energy that crosses the system boundary** is
  the orbit's: the Sun's input is a declared source in the module's own
  budget (`solar_input_across_system_boundary`).
- **So:** `fixed_surface_boundary_transfer` is a declared source or sink in
  the slab's own energy budget (`add_source` or `add_sink`, by its sign),
  named in its contract's sources and sinks.
- **What is unchanged:**
  - it is published per region;
  - the reconciler closes the system including it;
  - G2's requirement that it captures all the holding energy.
- **No reservoir module is added.**

Everything else in the protocol is unchanged.

## Erratum 2 to the code section, before any run (nothing above this line changed)

**Reading 1 (above).** `prescribed_temperature` goes on `SlabSurfaceModule`.

**Reading 2, from the configuration.** The tropical column's surface is
`ColumnSurfaceModule`, the slab beneath a radiating column; `SlabSurfaceModule`
is the airless slab. So the option is on `ColumnSurfaceModule`, under the same
name and with the same behaviour. It is passed through
`reference_column(prescribed_surface_temperature=...)` and
`two_column(prescribed_surface_temperatures=(rising, sinking))`. A held surface
starts at its held temperature, and a start elsewhere is refused.

**The tests are in `tests/integration/test_fixed_surface.py`:**

- the held temperature, exactly, every step;
- the holding contribution equal to minus the others;
- the published power, and the declared budget source, matching it;
- the reconciler passing;
- option-off identity;
- the refused start.

**The full suite passes: 676 tests.**
