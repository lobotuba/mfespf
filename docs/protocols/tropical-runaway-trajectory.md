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
