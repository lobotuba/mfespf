# Protocol: does N0's own trajectory run away for the ceiling's reason?

**Frozen before any run of this protocol.** Results are appended below a line
at the end; nothing above it changes afterwards. The design is the user's
(2026-10-14), recorded here with its implementation fixed.

**It is diagnostic only:**

- no model code changes;
- no defaults change;
- nothing is fixed.

## Prior results, not redone

**From `tropical-runaway-ceiling.md`, on prescribed states:**

- At a free-tropospheric RH of 0.95, the tropical column has no equilibrium
  inside RRTMG's tables. N = ASR − OLR stays above the band's 42.51 W m⁻²
  export, from 77.5 down to 52.1.
- At the observed humidity it equilibrates near 298.75 K.
- The stability threshold lies between r = 0.7 and 0.8 above 850 hPa.

**From `free-troposphere-moisture-budget.md`:** nothing in the single column
dries the free troposphere.

**The bulk-plume branch,** an attempt to supply that drying, is closed
(`convective-cape-jump-census.md`).

## Question

Does the prognostic N0 trajectory fail for the same radiative-moisture
reason? Or does another process contribute materially before that happens?

## Hypotheses

**Primary: N0 runs away because its free troposphere stays near RH 0.95.**
That makes the column's radiative requirement exceed the 42.51 W m⁻²
export throughout the valid temperature range.

**The alternatives, the only ones tested:**

| | alternative |
|---|---|
| **B, surface energy** | The surface or turbulent flux response, not the atmospheric humidity, supplies or fails to remove the destabilising energy |
| **C, adjustment heating** | Dry, convective or water adjustment contributes a material positive energy tendency, as net column energy and not vertical redistribution |
| **D, radiation domain** | The runaway begins only after the state leaves RRTMG's domain |

## The run

**N0, exactly as `ocean-heat-transport-circulation.md` ran it:**

- `oht_run.py n0`, the ±20° `tropical_column`;
- the declared physics, from the US Standard Atmosphere initial temperatures;
- master step 3600 s;
- up to 100000 steps.

**It continues as far as the existing refusal rules allow.** The earlier N0
stopped at step 45050, when the solver refused a negative optical depth.

**An analysis script reads the run's own retained trajectory.** It adds no
module and no output.

**G0, reproduction.** The rerun reproduces `oht_n0b.json` (the earlier N0's
tail-mean result, steps 36041-45050) to every printed digit. If it stops at a
different step, that is recorded and G0 fails.

## R1: the energy history, per step, then as daily means

**The state:**

- T_s and T(p);
- q_v(p) and RH(p), the published `atmosphere_relative_humidity`;
- **RH_300-700:** the mass-weighted mean RH over the layers centred between
  300 and 700 hPa;
- **RH_850-tp:** the mass-weighted mean RH over the layers centred between
  850 hPa and the cold point, the coldest layer below 50 hPa. That is the
  region the ceiling's r applied to.

**The top of the atmosphere:**

- ASR = `toa_downward_radiation` − `reflected_shortwave_radiation`;
- OLR;
- N = ASR − OLR;
- **N_res = N − 42.51.**

**The export, as applied:** `lateral_energy_export`, with
`lateral_latent_energy_export` reported separately.

**The surface:**

- net shortwave (`surface_absorbed_shortwave`);
- net longwave (`surface_net_radiative_flux` minus net shortwave);
- sensible and latent heat;
- the convective exchange with the air (`surface_convective_heat_flux`);
- the storage tendency, C_s·dT_s/dt from `surface_heat_capacity`.

**The atmosphere's adjustments,** column-integrated in W m⁻² from the
labelled temperature tendencies times layer heat capacity:

- `dry_convective_adjustment`;
- `latent_heat_release`;
- `surface_sensible_heat`;
- `tropical_exchange`;
- `radiative_flux_divergence`.

The vapour labels are reported as kg m⁻² day⁻¹:

- `convective_vapour_transport`;
- `vapour_eddy_diffusion`;
- `evaporation_minus_condensation`;
- `brewer_dobson_upwelling`.

**The radiation domain, per step:**

- the largest |`rrtmg_reference_temperature_departure`| over the layers,
  against RRTMG's ±30 K;
- `rrtmg_negative_optical_depths_floored`.

**The reconstruction, per step.** The column's storage rate is
d/dt[`atmosphere_column_energy` + C_s·T_s + L·`atmosphere_precipitable_water`],
with L the water module's reference latent heat. It is compared with
N − (`lateral_energy_export` + `lateral_latent_energy_export`), and the
residual reported in daily means.

**The run's conservation reconciliation is the formal energy gate.** It is
its worst relative residual over the completed steps. The reconstruction is
an independent check and is reported.

## R2: the chronology (from daily means)

| time | definition |
|---|---|
| **t_RH** | the first day from which every later daily mean of RH_300-700 is ≥ 0.80, the upper end of the ceiling's stability threshold |
| **t_N** | the first day from which every later daily mean of N_res is > 0, while the run is valid |
| **t_run** | the first day whose daily-mean T_s exceeds 301.61 K (the observed band's 299.61 K + 2 K), after which every 10-day mean of T_s exceeds the one before it, to the end of the valid data |
| **t_invalid** | the first step at which any layer's \|dT_ref\| exceeds 30 K, as a day |

**Also recorded:**

- the first step with a floored negative optical depth;
- the refusal step.

**Whether N_res ever reaches ≤ 0** in a daily mean before t_invalid, and the
lowest daily N_res reached. A stable column must reach it as it warms.

## Verdicts, frozen

**The existing radiative-moisture attribution is confirmed if, while the run
is valid (before t_invalid):**

1. RH_300-700 is ≥ 0.80 in its daily means from t_RH onward, and t_RH comes
   before t_run;
2. N_res never reaches ≤ 0 in a daily mean;
3. t_N comes no later than t_run;
4. B and C below are not supported.

| | supported only if |
|---|---|
| **B** | (i) the reconstruction fails: the daily-mean storage rate differs from N − export by more than 0.5 W m⁻² on more than 5 percent of valid days, and the discrepancy lies in the surface terms; or (ii) N_res ≤ 0 for 10 or more consecutive valid days while the 10-day-mean T_s keeps rising, so the surface warms by drawing on the atmosphere |
| **C** | the conservation reconciliation's worst residual exceeds 1e-12 over the valid steps, so some operator creates energy; or `dry_convective_adjustment` and `surface_convective_adjustment` together are a net column energy source of more than 0.5 W m⁻² in a daily mean. Otherwise each adjustment is reported as the redistribution it is. |
| **D** | t_invalid comes before t_N, or before the first daily N_res > 0 after spin-up. If t_invalid comes after t_run, the domain departure is recorded as a consequence, not the cause. |

**The attribution is overturned if** N_res reaches ≤ 0 and stays bounded
while the run is valid, and the runaway appears only through B, C or D. That
reopens the diagnosis.

**Anything else** is reported as found, with the times and values.

## What follows, frozen

- **If confirmed:** no isolated tropical-column correction is built. The
  next question is architectural: how ascent and subsidence can develop
  together from the start, so that the subsiding branch dries the free
  troposphere before either column reaches its moist runaway. It gets its
  own protocol, the user's decision.
- **If overturned:** the runaway diagnosis reopens.

## Recorded alongside: the global column (the user's wording)

> Free-tropospheric humidity is not considered internally closable in the
> global-mean single column. The experiments indicate that the missing
> physics is a drying pathway associated with spatial circulation, not
> another local humidity adjustment. No tuned free-tropospheric drying or
> humidity relaxation will be introduced solely to stabilize the single
> column.

---

## Erratum, before the run (nothing above this line changed)

Found in a 72-step smoke test of the analysis, before N0 was run. Both
readings are kept, and the verdicts are reported under both.

### 1. The export the column actually loses is about 45.5 W m⁻², not 42.51

**What the code does.** The tropical exchange's `lateral_energy_export` is
the prescribed band heat export, 42.51 W m⁻², plus the Brewer-Dobson
upwelling's net export of dry static energy (`tropical_exchange.py`). That is
45.53 W m⁻² in the smoke test. The latent counterpart,
`lateral_latent_energy_export`, is negligible there.

**The two readings:**

- **As frozen:** N_res = N − 42.51.
- **Reading 2, energy-consistent:** N_res,applied = N −
  (`lateral_energy_export` + `lateral_latent_energy_export`). The column
  equilibrates where this is zero.

t_N and verdict items 2 and 3 are evaluated under both.

**What this means for `tropical-runaway-ceiling.md`.** Its criterion N =
42.51 understated the column's sink by about 3 W m⁻². It does not change its
K1: H95's N stayed at 52.1 or above. Its HOBS crossing would move down by
about one 2.5 K step. It is recorded here and not re-run.

### 2. The frozen storage reconstruction does not close

**The formula.** d/dt[`atmosphere_column_energy` + C_s·T_s +
L·`atmosphere_precipitable_water`] misses N − export by 1.5-20 W m⁻² in the
smoke test's daily means. The likely causes are a missing store or the
outputs' differing time semantics. It would make B(i) fire falsely.

**The two readings:**

- **As frozen:** reported, as computed.
- **Reading 2:** the storage rate is the sum of the conservation reconciler's
  own interior energy stores, (stored_final − stored_initial)/dt for every
  energy budget except the system boundaries (orbit, space, lateral). In the
  smoke test it matches N − export to 2e-11 W m⁻².

B(i) is evaluated under both.

**Nothing else changes:** the run, the times, the thresholds and the verdict
rules.

## Result (appended; nothing above this line changed)

**The run.** N0 rerun exactly (`oht_run.py n0`).

- **It completed all 100000 steps.** Status: invalid on the validity-domain
  gate.
- **Worst relative residual:** 7.8e-16.
- **No refusal.** At step 45050, where the original N0 was refused, the
  current code floors and counts a negative optical depth
  (`rrtmg-domain-and-refusals.md`). The first floored step is 45050 exactly.

### G0: fails as frozen

The run did not stop at 45050, so its tail means (the last fifth of 100000
steps) are not the original's (steps 36041-45050). Of 199 shared printed
values, 139 differ.

- **Why it is likely the same trajectory to step 45050:** the first negative
  optical depth falls at exactly the original refusal step.
- **A supplementary check, not frozen:** a 45050-step rerun, whose tail is
  exactly the original's window, has been started. It is reported when it
  completes.
- **What the verdict uses:** only days 0-457 (steps up to about 10968), far
  before either.

### R2: the chronology (daily means)

| time | day | T_s | what happens |
|---|---:|---:|---|
| t_N, all three readings | **0** | 263.8 K, the cold start | R is positive from the start and stays positive through every valid day |
| t_RH | **197** | 298.8 K | RH 300-700 ≥ 0.80 on every later day. It first reaches 0.82 by day 50 and stays 0.81-0.93. |
| t_run | **308** | 301.7 K | T_s passes 301.61 K, and 10-day means rise thereafter |
| t_invalid | **457** | 303.9 K | the first \|dT_ref\| > 30 K |
| the first floored optical depth | step 45050, day 1877 | 315.2 K | |

**Order: t_N (0) ≤ t_RH (197) < t_run (308) < t_invalid (457).**

### The residual, three readings

| reading | sink | lowest daily value while valid | positive from |
|---|---|---:|---|
| **R_applied** (interprets) | N − applied export: 45.19-45.97 W m⁻², latent part ≤ 0.003 | +4.01 W m⁻² (day 1) | day 0 |
| R_45.5 | N − 45.5 | +4.05 | day 0 |
| R_frozen | N − 42.51 | +7.04 | day 0 |

**R_applied never reaches zero while the run is valid.** It falls as the
column warms but stays positive:

| day | T_s | RH 300-700 | ASR | OLR | R_applied |
|---:|---:|---:|---:|---:|---:|
| 10 | 282.8 K | 0.68 | 295.4 | 219.1 | 30.6 |
| 100 | 294.5 K | 0.83 | 301.2 | 234.3 | 20.9 |
| 197 | 298.8 K | 0.87 | 303.8 | 244.5 | 13.5 |
| 308 | 301.7 K | 0.88 | 305.6 | 250.1 | 9.9 |
| 457 | 303.9 K | 0.88 | 307.3 | 255.8 | 6.3 |

**After t_invalid it continues the same way, outside RRTMG's tables:**

| day | T_s | RH 300-700 | R_applied |
|---:|---:|---:|---:|
| 1000 | 309.4 K | 0.95 | 4.3 |
| 1877 | 315.2 K | 0.94 | 3.8 |
| 4165, the last | 320.1 K | 0.89 | **0.27** |

### B, C and D

**B: not supported under reading 2.**

- **The reconciler-store reconstruction** matches N − export to 6.1e-11
  W m⁻² on every valid day.
- **No 10-day period with R ≤ 0:** B(ii) does not fire.
- **The surface is near balance after spin-up.** Surface storage is about
  0-10 W m⁻². The latent flux rises from 68 to 160 W m⁻² as the column warms.

**The frozen reconstruction,** as the erratum anticipated:

- it misses by up to 20.5 W m⁻²;
- it exceeds 0.5 W m⁻² on 49 percent of valid days;
- so under the frozen reading B fires.

**C: not supported.**

- **Worst reconciliation residual:** 7.8e-16.
- **The dry adjustment and its surface counterpart together:** a net of at
  most 6.2e-11 W m⁻².
- **The dry adjustment's +37 W m⁻² into the air** (a mean over valid days)
  is the surface's convective exchange: redistribution, not creation.
- **The latent heat release, +143 W m⁻²,** converts the surface latent flux.

**D: not supported.** t_invalid (457) comes after t_N (0) and after t_run
(308). The departure from RRTMG's domain is a consequence of the warming, not
its cause.

### Verdicts

| reading | confirmed | B | C | D |
|---|---|---|---|---|
| as frozen | no | yes | no | no |
| erratum reading 2 | **yes** | no | no | no |

**As frozen,** the attribution is not confirmed, solely because B fires on
the reconstruction the erratum found defective before the run. Items 1-3
hold under the frozen R.

**Under the erratum's reading, the existing radiative-moisture attribution
is confirmed on N0's own trajectory:**

- **The humidity:** the free troposphere reaches the high-RH regime and stays
  in it, RH 300-700 ≥ 0.80 from day 197 and 0.81-0.93 from day 50.
- **The residual:** R_applied is positive from the start and never reaches
  zero while the radiation is valid.
- **The order:** the positive residual precedes the runaway (day 308), which
  precedes the domain departure (day 457).
- **No other source:** no surface, adjustment or domain process supplies
  the energy. The reconciler's stores account for all of it.

### What the result adds

**N0 is not an unbounded runaway.**

- **It warms slowly toward a hot equilibrium.** R_applied falls from 30 to
  0.27 W m⁻² by day 4165, at T_s = 320.1 K.
- **That equilibrium lies outside RRTMG's tables.**
- **It agrees with the ceiling's prescribed-state result.** There the H95
  column's N was 45.1 at an extrapolated 317.5 K, against the applied sink of
  about 45.5. That puts its crossing near 318-320 K, as N0's trajectory
  finds.

**The chain on the prognostic trajectory:** no free-tropospheric drying → RH
about 0.85-0.95 → OLR too small for the absorbed sunlight → R_applied > 0 →
warming toward an equilibrium near 320 K, beyond the radiation domain.

**What follows, as frozen:** no isolated tropical-column correction. The next
question is architectural, the user's decision: ascent and subsidence
developing together from the start.
