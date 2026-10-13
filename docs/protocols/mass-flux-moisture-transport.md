# Protocol: mass-flux convection with compensating subsidence

**Frozen before any code or run.** Results are appended below a line at the
end; nothing above it changes afterwards. The design is the user's
(2026-10-13 review). The parameter choices and their sources are declared
here, before any outcome is seen.

## Why

**`free-troposphere-moisture-budget.md`:** over 300-700 hPa,

- convection's vapour mixing supplies 96 percent of the water;
- condensation removes 98 percent;
- nothing in the column descends and dries;
- denying the supply made the free troposphere wetter: 0.985, Db.

**`tropical-runaway-ceiling.md`:** the radiation is right, and the free
troposphere must stay below about 0.75 for the tropics to be stable.
R4 holds 0.89.

**What the code does.** Convection is a lapse-rate adjustment that mixes
specific humidity uniformly over each mixed group. Real convection rises in
narrow plumes, and the rest of the column sinks to compensate. The sinking
air dries, and only detrainment and evaporating condensate re-moisten it
(Yanai et al. 1973, *J. Atmos. Sci.* 30, 611; Arakawa & Schubert 1974,
*J. Atmos. Sci.* 31, 674).

## The hierarchy

dilute CAPE → cloud-base mass flux M_b → plume M_u(z) → M_e(z) = -M_u(z) →
conservative scalar fluxes

A new convection option, `mass_flux_convection`, replaces the moist,
entraining-plume lapse-rate adjustment and its vapour mixing. The dry
adjustment of dry-unstable layers and the surface sensible heat flux are
unchanged.

### 1. Trigger and amplitude, separately

- **The source layer:** the layer of largest moist static energy
  h = c_p T + g z + L q among those centred below 850 hPa. That is ZM's
  choice of launch level (Zhang & McFarlane 1995, *Atmos.-Ocean* 33, 407).
- **Dilute CAPE:** CAPE_d is the integrated buoyancy of the entraining plume
  of section 2, launched from the source layer, from its level of free
  convection to its termination level.
- **The trigger:** convection runs only if CAPE_d > C_min = 70 J kg^-1.
  That is CAM5 `zm_conv.F90` `capelmt`.
- **The closure:**

  M_b = max[0, (CAPE_d - C_residual) / (τ_c Γ_M)]

  - **C_residual = 0,** as ZM removes CAPE.
  - **τ_c = 3600 s,** CAM5 `zm_conv.F90` `tau`. It is a declared closure
    hypothesis, not a numerical constant (section 7).
  - **Γ_M:** the rate at which unit cloud-base mass flux removes CAPE_d,
    diagnosed with the model's own operator. Apply the transport operator
    of section 4 at M_b = 1 kg m^-2 s^-1 for one step dt, recompute CAPE_d
    on the result, and Γ_M = -(ΔCAPE_d)/dt.
  - **If Γ_M ≤ 0,** M_b = 0 and the step is counted.

### 2. The plume

- **Continuity:** ∂M_u/∂z = (ε - δ) M_u, with M_u(z_b) = M_b at the level
  of free convection z_b. Below z_b, M_u ramps linearly in pressure from
  zero at the source layer's base.
- **The entrainment law: ε constant, 1 km^-1.**
  - The source is Romps (2010, *J. Atmos. Sci.* 67, 1908), direct LES
    measurement in deep convection.
  - The declared range is 0.5-2 km^-1: Romps finds bulk methods
    underestimate entrainment by about a factor of 2.
  - The same ε dilutes the CAPE parcel. Neale et al. (2008, *J. Climate* 21,
    5904) use 0.69 km^-1, inside the range.
- **The detrainment law:** δ = ε below the termination level. That is the
  turbulent exchange, with constant mass flux (Tiedtke 1989, *Mon. Wea. Rev.*
  117, 1779). The whole remaining M_u detrains in the termination layer, as
  organised detrainment.
- **Plume thermodynamics:**
  - h_u and q_t,u are integrated upward with entrainment of environmental
    h_e and q_e;
  - the plume is saturated above its lifting condensation level, at its
    own temperature;
  - condensate formed in a layer leaves the plume at that layer's top, as
    P_f.
- **The buoyancy:**

  B_u = g[(T_u - T_e)/T_e + 0.61(q_v,u - q_v,e) - q_c,u]

  q_c,u is the condensate formed within the layer before it leaves: the
  loading, at this resolution.
- **Termination (simplification, declared):** the plume stops at the level
  where its integrated buoyancy work from z_b, ∫B_u dz, first returns to
  zero or below after its level of neutral buoyancy. A kinetic-energy
  equation for w_u is deferred.
- **Downdrafts are deferred:** M_d = 0.

### 3. Compensation by exact continuity

- **The environment's flux:** at every interface, M_e = -(M_u + M_d) =
  -M_u. Subsidence has no amplitude of its own.
- **The rule (A1):** a step whose interface mass continuity misses by more
  than round-off is rejected.

### 4. One conservative operator

- **The flux.** For each conserved scalar φ, moist static energy h and total
  water q_t:

  F_φ = M_u (φ_u - φ_e)

- **The update:** (∂φ̄/∂t)_conv = -g ∂F_φ/∂p + S_φ.
  - **S_φ holds only** the precipitation that leaves, and the evaporation of
    section 6.
  - **Entrainment, detrainment, plume transport and subsidence** are this one
    operator. None is applied separately.
- **Numerics (declared).** Fluxes are upwind at interfaces. If M dt exceeds
  half of any layer's mass, the operator subcycles until it does not. M_b is
  never limited: CAM5's Courant cap on M_b would be a numerical limit on a
  physical quantity, which CLAUDE.md forbids.

### 5. The non-targeting rule

- **What may not set M_b or subsidence:** neither may be closed on, or
  tuned against, relative humidity, precipitable water, precipitation or
  evaporation.
- **What may not move once frozen:** τ_c, ε, the detrainment law, α, f_evap
  and C_min. Their sensitivity is preregistered in section 7.
- **What stays a prediction:** RH, precipitable water and precipitation.

### 6. α = 0.67

- **The measured `EVAPORATING_CONDENSATE_FRACTION`** partitions what forms:
  - P_evap = α P_formed;
  - P_surface = (1 - α) P_formed.
- **The deposition rule (frozen):**
  - P_evap evaporates into the environment in the layer where it formed,
    up to that layer's saturation deficit. That is Romps' (2014)
    placement.
  - What exceeds the deficit falls to the next layer below and evaporates
    there, again limited by subsaturation, and so on downward.
  - What reaches the surface joins P_surface.
- **The energy:** the latent cooling is the environment's.
- **α does not close M_b.**

### 7. τ_c as a hypothesis, and the preregistered sensitivity

**The verdict runs use τ_c = 3600 s and ε = 1 km^-1.**

**The sensitivity ensemble,** the same runs at:

- τ_c = 1800, 7200 and 14400 s. Observed dilute CAPE implies a relaxation
  of 3-5 h, per the ARM-ASR literature cited in the E3SM ZM documentation;
- ε = 0.5 and 2 km^-1.

**A qualitative conclusion holds only if it survives the whole ensemble.**
Nothing is retuned after the results.

**The alternative, named and excluded.** A closure on the subcloud convective
velocity scale, M_b = C_M ρ_b w_*, is a named alternative hypothesis for
shallow regimes. It is not part of this protocol.

## Ownership and bookkeeping

**Who owns what.** Convection computes the operator and publishes the
temperature and vapour tendencies, the precipitation reaching the surface,
and the evaporation. The water module remains the owner of vapour and of
latent energy:

- it commits the vapour tendency;
- it records the precipitation transfer to the surface water store;
- it books the latent energy of the plume's condensation, which leaves
  (1 - α) as precipitation;
- it books the latent energy of the evaporation, as transfers both modules
  declare.

**What stays.** The water module's surface evaporation, eddy diffusion and
end-of-step saturation are unchanged and compose after convection, as now.

## The gates (all must pass before any physics run)

| gate | passes if |
|---|---|
| A1: mass | at every interface M_u + M_e = 0 to round-off; entrainment and detrainment exchanges cancel |
| A2: water | per step, column water change = E - P_surface to 1e-12 relative; plume, environment and evaporation transfers contribute zero |
| A3: energy | column moist enthalpy closes to the existing reconciliation tolerance, phase changes and precipitation enthalpy included |
| A4: drying without destruction | a moist, conditionally unstable column with surface fluxes off: total water and moist enthalpy constant less what precipitates; the subsiding environment warms and its RH falls; no water leaves except as precipitation |
| A5: no circulation | with CAPE_d ≤ C_min, M_u = M_e = 0 and the operator changes nothing. With the option off the column is the present one, to the bit (R4, 200 steps) |
| A6: timestep | from R4's state, 10 days at dt = 3600, 1800 and 900 s: free-tropospheric RH, precipitation, M_b and the energy residual converge (900-vs-1800 differences at most half the 1800-vs-3600 ones) |
| A7: Db reversed | from R4's state, halving the surface evaporation coefficient lowers the free troposphere's RH: less supply now means drier, where Db's denial made it wetter |

## The tests of the physics

**P1, its own theory, judged.** In layers between 300 and 700 hPa where
M_u > 0 and δ = ε, the equilibrated RH matches Romps' (2014) environmental
budget with the run's own values, within 0.05:

RH = (a + α)/(1 + a), a = ε(1 - α)/γ

- γ is the run's water-vapour lapse rate, `water_vapour_lapse_rate`;
- the in-repo implementation is `bulk_plume_relative_humidity` with δ = ε;
- the scheme must produce the humidity its own physics predicts.

**P2, the global column:** R4's configuration with `mass_flux_convection`,
200000 steps.

- **Judged:** G0-G4 as frozen in `global-column-chou-longwave.md`.
- **Reported against R4:**
  - RH 300-700;
  - precipitable water;
  - the C-shape;
  - the lapse rate;
  - cloud-top temperatures;
  - M_b and its distribution;
  - the trigger frequency;
  - the CAPE_d distribution.

**P3, the tropical column:** N0's configuration
(`ocean-heat-transport-circulation.md`) with `mass_flux_convection`.

- **Judged:** valid, and equilibrated (net TOA within ±0.1 W m^-2 of its
  export).
- **Reported:** the surface, against 299.61 K.

**The sensitivity ensemble of section 7** runs on P2 and P3.

## What is expected, stated before any code

**The theory has a floor.** With α = 0.67, Romps' RH = (a + α)/(1 + a) can
never fall below α: a single convecting column under this theory holds its
free troposphere at 0.67 or above, whatever the entrainment.

- **At ε = 1 km^-1 and γ about 0.45 km^-1:** a ≈ 0.73 and RH ≈ 0.81.
- **Expected outcome:**
  - the first subsaturated free troposphere the column has had, and drier
    than R4's 0.89;
  - above the 0.75 tropical threshold, so P3 may still run away;
  - G4 expected to fail.
- **The observed 0.43 is out of reach for a single column.** Reaching it
  needs the dry subtropics: the two-column overturning, on a base state that
  no longer runs away. That is the step after this one.

**Why build it anyway.** It replaces the closure the budget shows to be
structurally unable to dry. It makes subsidence a consequence of convection
rather than a parameter. It is the convection every later configuration
needs. P1 makes it falsifiable against its own theory.

## Order

1. The code, with A1-A5 as unit tests;
2. A6 and A7;
3. P2 and P3;
4. P1 on both;
5. the sensitivity ensemble.

**If any of A1-A5 fails,** the build stops and the failure is recorded.

**Adoption is the user's decision.**

---

## Erratum, before any code (nothing above this line changed)

Mapping the frozen design onto the code found three places where the text
cannot be implemented as written, or is ambiguous. Each is recorded with both
readings. No outcome has been computed.

### 1. Γ_M at unit mass flux for a whole step

**What was frozen:** "apply the operator at M_b = 1 kg m^-2 s^-1 for one
step dt". At dt = 3600 s that moves 3600 kg m^-2: more than the whole column
above 500 hPa. The CAPE change it gives is that of a column turned over many
times, not a rate.

**The readings:**

1. **As frozen,** which is not a rate, so it is not taken.
2. **Γ_M as the derivative it names,** taken. It is the rate at which unit
   cloud-base mass flux removes CAPE_d.
   - The operator's tendencies are linear in M_b. They are formed at
     M_b = 1, as rates.
   - CAPE_d is re-evaluated on the state advanced along them by
     M_b Δt = 1 kg m^-2: a small displacement, under 1 percent of any
     tropospheric layer's mass.
   - Γ_M = -ΔCAPE_d / (1 kg m^-2 / (1 kg m^-2 s^-1)) is the finite
     difference of that derivative. The trigger, the closure form, τ_c and
     C_residual are unchanged.

### 2. Who heats the air with the plume's condensation

**The ownership.** The water module owns vapour and latent energy, and it is
where every phase change's latent heat enters (CLAUDE.md, "a process that
moves matter must move its energy in the same interval").

**The operator's net vapour change** in layer k splits exactly into:

- a conservative transport t_k (column sum zero);
- the net convective condensation c_k = P_f,k - e_k (column sum
  P_surface).

**The readings:**

1. **Convection applies the latent heating itself.** That moves a phase
   change's energy out of the module that owns it.
2. **Convection publishes t_k and c_k,** taken.
   - The water module commits the vapour change t_k - c_k.
   - It adds c_k to its condensation: precipitation, latent heating L c_k,
     and the fusion excess where c_k is deposition. That is the same route
     its own condensation takes.
   - It includes L c_k / C_k in the end-of-step temperature its saturation
     adjustment uses.
   - Convection's own temperature tendency then excludes L c_k. It
     redistributes the column's enthalpy and sums to zero, as the adjustment
     it replaces did.

### 3. The dry adjustment's vapour mixing

**The ambiguity.** The frozen text replaces "the moist, entraining-plume
lapse-rate adjustment and its vapour mixing" and keeps "the dry adjustment
of dry-unstable layers". It does not say whether the dry groups still mix
vapour.

**The readings:**

1. **They do not,** which leaves the boundary layer's turbulence moving heat
   and not water.
2. **They do,** taken. Mixing within a dry-unstable group is turbulence, and
   it carries water as it carries heat.
   - With the dry-adiabat target, dry-unstable groups sit near the surface.
     The free troposphere's water then comes only through the plume and its
     subsidence.
   - That vapour mixing is reported separately, so its share stays visible.

### Numerics, made precise

The subcycling criterion counts every outflow from a layer in a substep, not
only M:

- the subsiding environmental flux out through its base;
- entrainment into the plume.

Each substep keeps their sum, times the substep, at no more than half the
layer's mass. M_b is still never limited.

## Second erratum, before any code (nothing above this line changed)

**What was frozen:** M_u = M_b at the level of free convection z_b, ramping
linearly in pressure from zero at the source layer's base.

**Why it cannot stand.** The plume's composition below z_b would then depend
on the ramp's shape, which depends on where z_b is. z_b depends on the
plume's buoyancy, which depends on its composition. That is circular. The
exact discrete budgets also need one plume that both sets CAPE_d (section 1)
and carries the transport (section 4). A ramped plume and the dilute CAPE
parcel are two different plumes.

**The readings:**

1. **The ramp, iterated to a fixed point over z_b.** That needs a further
   numerical rule nobody declared, and may not converge.
2. **Cloud base at the launch level, as ZM has it,** taken (Zhang & McFarlane
   1995).
   - The plume takes all its mass, M_b, from the source layer: E = M_b and
     D = 0 there.
   - From the source layer's top it rises at constant M_b, with δ = ε above.
   - The whole mass detrains in the termination layer.
   - The dilute CAPE parcel and the transported plume are then the same
     plume, layer for layer. CAPE_d is integrated from z_b, the first
     positively buoyant layer above the source, as frozen.

**The plume's discrete form, exact for conservation:**

- **Within a layer,** the plume mixes its inflow with entrained air,
  conserving moist static energy and total water exactly.
- **Above saturation,** it then condenses at the layer's pressure:
  - h is unchanged;
  - the condensate formed is P_k;
  - it leaves the plume at the layer's top.
- **The environment's operator,** per unit M_b:

  m_k dφ_k/dt = M_(k+1) φ_(k+1) - M_k φ_k - E_k φ_k + D_k φ_u,k

  for φ = s (dry static energy) and q. The first two terms are the
  compensating subsidence, upwind.
- **Exactness.** Entrainment draws the start-of-step environment, so the
  plume's budget identities, Σ(E h_e - D h_u) = 0 and Σ(E q_e - D q_u) = ΣP,
  hold exactly under subcycling.

## Third erratum, before any gate result (nothing above this line changed)

A development probe ran the scheme from R4's state. It is not a gate, and its
numbers decide nothing. It found two defects in the numerics as frozen.

1. **The closure was applied once per step.** M_b = CAPE_d / (τ_c Γ_M) is a
   rate law. Held for the whole 3600 s step, it extrapolates linearly. Where
   Γ_M is small, that is a mass flux that overturns the column many times in
   one step, long after its CAPE is gone.
2. **Entrainment drew the start-of-step environment in every substep** (the
   second erratum). Under many substeps it kept removing water at the start
   value from a layer that subsidence had already refilled with drier air.
   The probe refused at step 54: "mass-flux transport drew a layer's vapour
   below zero".

**The readings:**

1. **Cap M_b,** as CAM5 does with its Courant limit. That is a numerical
   limit on a physical quantity, which the protocol forbids, so it is not
   taken.
2. **Integrate the closure as the rate law it is,** taken.
   - Each substep recomputes the plume and the closure (CAPE_d, Γ_M, M_b) on
     the current state.
   - Each substep then takes one exact, conservative step of the operator,
     long enough that no layer loses more than half its mass.
   - Entrainment draws the current environment, so the plume's budget
     identities hold exactly within every substep.
   - The measured fraction of that substep's condensate evaporates on the
     current state's deficits.
   - Convection stops for the rest of the step once CAPE_d falls to the
     trigger, or Γ_M to zero or below.
   - A step that would need more than 500 substeps is refused, naming the
     state, rather than limited.

The trigger, τ_c, ε, α, the closure form and every gate are unchanged.

## A1-A7 result (appended; nothing above this line changed)

**A1-A5 pass,** as unit tests in `tests/unit/test_mass_flux.py`, commit
4ef5e1f:

- **The plume's budgets:** dM = E - D; Σ(E h_e - D h_u) = 0 to 1e-12;
  Σ(E q_e - D q_u) = ΣP to 1e-12.
- **The operator:** the energy redistribution sums to zero to 1e-10 of a
  layer's change; the transport sums to zero; the net condensation sums to
  the surface precipitation; the column's water falls by exactly that.
- **Drying without destruction:** an hour of the closure on a moist,
  unstable column warms and dries the environment between the level of free
  convection and the detrainment layer, and water leaves only as
  precipitation.
- **No CAPE:** no convection. The evaporation cascade passes what a
  saturated layer cannot take to the layer below.
- **A coupled convecting column** reconciles its water and energy.
- **With the option off,** 200 steps of R4 give all 13000 field values
  identical to the code before the build, to the bit.

**Found and corrected while building, before any gate was judged:** the
water module subtracted convection's condensation twice, once from the
provisional profile and once as its own. The conservation gate refused the
first coupled run, and the cause is recorded in the commit.

**A6 and A7:** from R4's state, 10 days each. All four runs are invalid on
RRTMG's domain. The top layer, 1.1 hPa, drifts to 36 K below RRTMG's
reference; it starts at -27.7 K in R4's state. The water and energy
residuals are 5-6e-16.

| run | dt | RH 300-700, 10-day mean | precipitation, mm/day | mean M_b, kg m^-2 s^-1 | share of steps convecting | surface at day 10 |
|---|---:|---:|---:|---:|---:|---:|
| A6 | 3600 s | 0.9191 | 3.652 | 0.00635 | 0.34 | 309.63 K |
| A6 | 1800 s | 0.9170 | 3.678 | 0.00664 | 0.18 | 309.54 K |
| A6 | 900 s | 0.9204 | 3.668 | 0.00664 | 0.15 | 309.56 K |
| A7 | 3600 s, half the evaporation coefficient | **0.9313** | 3.409 | 0.00550 | 0.30 | 311.60 K |

**A6 fails as frozen.** The criterion is that the 900-vs-1800 s difference be
at most half the 1800-vs-3600 s one:

- precipitation passes (0.010 against 0.026);
- M_b passes (0.000 against 0.0003);
- **RH fails** (0.0034 against 0.0021): the differences are small but do not
  shrink;
- the energy residual is at round-off in all three (5.1, 5.2 and
  5.8e-16), so its "convergence" is the noise of the arithmetic. Read
  literally, it fails too.

**A7 fails as frozen.** Halving the surface evaporation made the free
troposphere *wetter*, 0.931 against 0.919, not drier. With less evaporative
cooling, the surface warms 2 K more in the 10 days. This is a drifting
transient: R4's state warms about 8 K in 10 days under this scheme. Db's
pathology has not been reversed.

**As frozen, the gates must all pass before any physics run, so P1-P3 and the
sensitivity ensemble do not run.** The build stands, conservative and
tested, and switched off by default.

**Reported, for the user's decision.** Across the probe and these runs:

- **The scheme convects intermittently,** at the realistic magnitude of a
  ZM-type scheme: M_b of 0.01-0.03 kg m^-2 s^-1 when active, about 0.006 on
  average, which is 54 hPa/day of mean subsidence.
- **The free troposphere stays near RH 0.92.** That is not the 0.81 Romps'
  relation predicts for its ε and α.
- **The likely reason,** not measured:
  - with ε = 1 km^-1, dilute CAPE above the 70 J/kg trigger exists only in a
    humid free troposphere (RH 0.8-0.95 on the synthetic soundings);
  - convection therefore runs only while the column is humid, and between
    events nothing sinks while radiation cools and moistens it.

  The trigger, not the plume's budget, then sets the humidity.
