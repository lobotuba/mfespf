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

## S1 result (appended; nothing above this line changed)

**Runs, as frozen.** 20000 steps did not meet the equilibrium criterion, so
the run was extended in 20000-step increments: 40000, then 60000, the
maximum. The verdict is read from the 60000-step run, over its last fifth.

**The equilibrium criterion is narrowly unmet even at 60000.** The
atmosphere's stored-energy rate is +0.104 W m⁻² in the rising column, just
outside ±0.1, and −0.023 in the sinking one. The descent RH's ten 50-day means
lie within 0.187-0.206 (range 0.019, within 0.02).

**The state barely moved over the extensions:**

| | 20000 steps | 40000 | 60000 |
|---|---:|---:|---:|
| descent omega500, hPa/day | 26.37 | 26.34 | 26.37 |
| descent RH 300-700 | 0.187 | 0.186 | 0.185 |
| band RH 300-700 | 0.607 | 0.606 | 0.606 |
| ΔF, W m⁻² | 142.7 | 143.1 | 142.7 |

### Gates

| gate | result | detail |
|---|---|---|
| G1 conservation | **pass** | reconciler residual 1.0e-15; every report passes |
| G2 the fixed temperature | **pass** | both surfaces exactly at 300.26 and 298.79 K at every step (deviation 0.0) |
| G3 the WTG exchange | **pass** | area-weighted cancellation 6.7e-16 relative (dry static energy) and 1.6e-16 (vapour) |
| G4 validity | pass | the largest \|dT_ref\| is 22.7 K, and the validity gate passes |
| G4 boundedness (1 percent drift between 50-day means) | **fail** | descent precipitable water 2.9 percent, descent dry adjustment 7.7 percent, descent latent heat release 349 percent (intermittent), ascent latent heat release 1.5 percent. All other quantities are under 1 percent. |
| G5 orientation | **pass** | descent omega500 positive on 100 percent of late days |

### G6, the descent column's free-tropospheric vapour budget (kg m⁻² day⁻¹, 300-700 hPa)

| label | late mean |
|---|---:|
| circulation (subsidence) | **−0.207** |
| convective vapour transport | +0.188 |
| eddy diffusion | +0.021 |
| evaporation − condensation | −4e-17, no condensation |
| Brewer-Dobson upwelling | 0 |

**The circulation is the only drying term,** and it balances the convective
and diffusive moistening.

### G7, the frozen predictions

| prediction | model | observed | passes if | |
|---|---:|---:|---|---|
| descent omega500 | 26.4 hPa/day | 14.3 | 7-29 | **pass** |
| descent RH 300-700 | **0.185** | 0.35 | 0.25-0.45 | fail, too dry |
| band RH 300-700 | **0.606** | 0.43 | 0.33-0.53 | fail, too wet |
| C-shaped descent profile | minimum 0.097 at 300 hPa | 0.28 at 400-500 | minimum within 350-650 hPa | fail |
| band precipitable water | 68.5 kg m⁻² | 38.0 | 28.5-47.5 | fail |
| T500, rising − sinking | +0.25 K | +0.65 | +0.33 to +1.30 | fail |

**The profiles at NCEP's levels:**

| level, hPa | 1000 | 925 | 850 | 700 | 600 | 500 | 400 | 300 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| descent RH, model | 0.82 | 0.87 | 0.95 | 0.23 | 0.29 | 0.19 | 0.11 | 0.10 |
| descent RH, observed | 0.75 | 0.69 | 0.57 | 0.39 | 0.32 | 0.28 | 0.28 | 0.50 |

**The rising column's RH 300-700 is 0.94,** against 0.50 observed.

### Reported, not gated

**The implied surface-transfer contrast:**

- ΔF = F_hold,rising − F_hold,sinking = **+142.7 W m⁻²** (+61.3 into the
  rising surface, −81.4 out of the sinking one);
- **about 10 times the observed ocean delivery contrast** (OAFlux 13.39, NCEP
  4.67);
- the sign agrees, the ocean delivering into the rising region.

**The band mean of F_hold is −1.76 W m⁻²:** a small net, so the holding is
almost all redistribution.

**Each column's TOA flux against its export:**

| | N | applied export |
|---|---:|---:|
| rising | 60.2 W m⁻² | 45.6 |
| sinking | 31.3 W m⁻² | 45.6 |

**The surfaces:**

| | latent heat flux | clear-sky OLR | all-sky OLR | cloud LW effect |
|---|---:|---:|---:|---:|
| rising | 159.9 W m⁻² | 263.7 | 244.5 | 19.2 |
| sinking | 93.5 W m⁻² | 298.0 | 269.1 | 29.0 |

The sinking column's evaporation is 58 percent of the rising column's.

### Verdict, as frozen: anything else, reported as found

- **Not "drying pathway supported":** descent RH and band RH are outside
  their ranges, and G4's boundedness fails.
- **Not "local moistening blocker":** the descent column is dry, not ≥ 0.75.
- **Not "surface contrast insufficient":** G5 holds and the subsidence is
  26.4 hPa/day.
- **Not "energetically unrealistic":** that verdict requires drying
  supported. Its ΔF condition (above 26.8 W m⁻²) would be met.

**S2 is therefore not opened,** as frozen.

### What it says, and what it does not

**1. The subsidence-drying pathway exists in this coupling.**

- Under the observed 1.47 K contrast, the relaxed spectral WTG circulation
  forms with the right orientation and realistic strength (26 hPa/day).
- The circulation alone dries the descent column's free troposphere, against
  convective moistening.
- It makes a moist boundary layer (RH 0.82-0.97 up to 828 hPa) capped by a
  sharp drop to 0.24 by 726 hPa: the shape of a trade-wind inversion.
- That is the first dry subsiding branch this model has produced.

**2. The partition is wrong in both columns.**

- **The rising column saturates its free troposphere (0.94).** That is the
  same convective vapour mixing `free-troposphere-moisture-budget.md`
  measured.
- **The sinking column overshoots dry (0.19 against 0.35).**
- **Its moist layer is capped near 830 hPa and saturated at the cap (0.97),**
  where the observed 850 hPa RH is 0.57. Nothing ventilates the boundary
  layer upward the way trade cumulus do.
- **Read with care:** the observed rising region's 0.50 mixes convecting and
  non-convecting air at 2.5° and monthly means
  (`ocean-heat-transport-circulation.md`). The rising column's comparison is
  not like for like, though the band's is.

**3. The observed contrast is not energetically sustainable here.**

- Holding it needs a surface-energy redistribution 10 times the observed
  ocean transport.
- The sinking surface gains energy: its evaporation is low under a capped,
  near-saturated boundary layer.
- The rising surface loses it.
- With prognostic surfaces under the observed 13.39 W m⁻², the contrast would
  be expected to collapse. That is a prediction, and S2 was not run.

**4. A declared limitation that bears on point 3.** Both columns carry the
same prescribed ISCCP tropical-mean cloud field. In reality the subsiding
regions' low clouds reflect sunlight and cool their surface. This protocol
cannot represent that.

**Not concluded:**

- no adoption;
- no default changed;
- no S1 state taken as a coupled solution.
