# Protocol: the fall-path closure candidate F0, under frozen radiation

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards.

**Authorised by the user (2026-10-13)** after
`convective-condensate-ledger.md`.

**The sequence:**

1. this protocol (F0 under frozen radiation, with a timestep study);
2. the stratospheric protocol;
3. F0 with interactive radiation;
4. reconsidering P1-P3.

**Constraints:**

- P1-P3 stay blocked.
- No α is adopted, imposed or used as an acceptance target.
- c0 and ke are not fitted. They are a baseline hypothesis.
- Passing here makes F0 "a numerically acceptable closure candidate awaiting
  radiative-coupling validation", not adopted physics.

## F0

Aα3 at the nominal constants, with every pathway physical and local, and no
α target:

- **cloud water** carried by the plume, converted implicitly at
  c0 = 0.045 m⁻¹ (CAM5 ZM `cldprp`, `c0_ocn` in its WRF port);
- **detrained cloud water** evaporating to the receiving layer's wet-bulb
  capacity. The rest becomes precipitation in that layer;
- **precipitation** falling interface by interface:
  - it evaporates under CAM5 ZM's rate law, ke = 1e-6;
  - each layer is limited by its remaining wet-bulb capacity and by what
    falls into it;
  - rain formed in a layer is not evaporated there;
  - the rest reaches the surface;
- **no redistribution** of evaporation a layer cannot take.

**The options:**

- `condensate_pathway="carried"`;
- `cloud_water_conversion_rate=0.045`;
- `precipitation_evaporation="rate_law_only"`;
- `detrained_state="plume"`.

**Everything else as Aα3:**

- R4's state and configuration with mass-flux convection;
- the common prescribed radiation, replayed;
- ε = δ = 1 km⁻¹, trigger 70 J/kg, τ_c = 3600 s.

**Stage A measured F0 at 3600 s** at α_total = 0.077. That is a recorded
result, not a gate.

## Molecule-once quantities: a provenance tracer

The user's design, frozen. Gross Eulerian phase-change totals (the existing
ledger) are kept separate from these.

### The history classes

Water carries an irreversible history:

U → C → R → (surface: the cohort ends)

- **U:** never condensed.
- **C:** condensed at least once, never formed precipitation.
- **R:** formed precipitation at least once.

**Physical phase changes within a class,** and its history does not. A
molecule that condenses, evaporates and condenses again enters C on its first
condensation and stays C.

**The cohort boundary is the surface.**

- Water that reaches the surface completes its cohort.
- Water the surface later evaporates starts a new cohort as U.
- Ocean-atmosphere cycles are not remembered.

**Provenance bins.** C and R carry the layer of first condensation, i:

- C_i and R_i, where i is a convective birth layer;
- U is one class.

**Large-scale condensation in the water module** rains out immediately, with
no evaporation:

- water it condenses completes the whole sequence U → C → R → surface in that
  step;
- it needs no persistent class;
- it is counted in its own accumulators, labelled by process and layer.

**The tracer state.** For each layer, the environmental vapour split into
U, C_i (i = 1..n) and R_i: 2n + 1 classes. The sum of the classes is exactly
the layer's vapour.

**Within a substep** (steady, as the plume is), the plume carries its vapour
and its cloud water by class. Precipitation falls by class. Nothing is stored
aloft between steps, so no stored cloud or rain classes are needed.

### Mixing and transport

Every operator applies to each class exactly as it does to the vapour:

- **subsidence, entrainment and detrainment:** the flux-form operator, per
  class;
- **the plume:**
  - its vapour and cloud water mix linearly with entrained air by class;
  - it detrains each class in proportion;
- **the dry adjustment's mixing:** uniform specific humidity of each class
  within each group, from the same base profile;
- **eddy diffusion:** the same linear operator, per class.

### Phase changes and milestones

Removal from a pool is in proportion to each class's share of that pool. That
is the declared well-mixed assumption.

| process | classes | milestone |
|---|---|---|
| condensation in the plume at layer k | vapour U → cloud C_k; C_i → C_i; R_i → R_i | **first_condensed[k] += U share** |
| evaporation of carried cloud in the plume | cloud → vapour, same class | none |
| conversion of cloud to precipitation at layer k | C_i → rain R_i; R_i → R_i | **first_precip_origin[i] += C_i share**, and **first_precip_location[k]** |
| detrained cloud water | evaporates into the receiving layer's vapour, same class; or becomes rain (C_i → R_i, a first precipitation at k) | first_precip, where it becomes rain |
| rain evaporation | rain R_i → vapour R_i in that layer | none |
| rain reaching the surface | R_i leaves | **first_surface_origin[i]** |
| large-scale condensation at layer k | removes each vapour class in proportion | U share: **first_condensed_ls[k]**, which is also its first precipitation and first surface arrival. C_i: first_precip_origin[i] and first_surface_origin[i]. R_i: first_surface_origin[i] |
| surface evaporation | enters U at the lowest layer | none |

**A layer the water module empties exactly** empties every class.

### Quantities, convective-born

All are over the reading window, with i the layer of first condensation:

- C_unique = Σ first_condensed;
- C_precip,unique = Σ first_precip;
- P_surface,unique = Σ first_surface;
- FE = C_precip,unique / C_unique;
- SE = P_surface,unique / C_precip,unique;
- PE = P_surface,unique / C_unique = FE·SE;
- α_total,unique = 1 − PE.

**Conditioned on origin:**

- FE(i) = first_precip_origin[i] / first_condensed[i];
- SE(i) and PE(i) likewise.

**Also reported:**

- the formation-location distribution, first_precip_location;
- large-scale-born and all-born water, separately;
- the inventory change of convective-born water over the window. The ratios
  are window fluxes, and that change says how far from steady they are.

**Internal names are milestones:** `first_condensed`, `first_precip`,
`first_surface`. They are never `C_precip`, so gross conversion cannot be
fed into them.

## The gates (frozen radiation)

**F1, conservation:**

- the conservation gate passes;
- the worst relative residual is ≤ 1e-12.

**F2, no negative states:**

- every run completes every step, with no refusal;
- no tracer class falls below −1e-12 of its layer's vapour.

**F3, no forced recycling.** R_recycle = C_recycle / (convective evaporation
into the environment), over the window.

- C_recycle is the water module's counterfactual measure of this step's
  recondensation of convective evaporation.
- **The condition:** R_recycle ≤ one tenth of E0's 0.042, so ≤ 0.0042, at
  every timestep.

**F4, timestep convergence.**

- **The runs:** dt = 3600, 1800 and 900 s, 20 days each (480, 960 and 1920
  steps).
- **The quantities:**
  - RH 300-700;
  - the RH profile, read as the largest difference over the layers from
    100 to 850 hPa;
  - surface precipitation;
  - FE, SE and α_total,unique;
  - cloud-base mass flux;
  - duty cycle;
  - CAPE_d;
  - subsidence per day;
  - R_recycle.
- **A6's frozen rule, unchanged:** the 900-vs-1800 s difference is at most
  half the 1800-vs-3600 s one.
  - It is applied to the late-period means (days 10-20), which decide.
  - It is applied to the full-run means and to the last day's means, which
    are reported only.
- **Each quantity is classified on the late means:**
  - **converged:** the rule holds;
  - **phase-sensitive:** the rule fails, but every pairwise difference is
    within twice the combined standard error of the two means. The standard
    error is the spread of daily means over √(number of days).
  - **changed equilibrium:** the rule fails and some difference exceeds
    that.
- **F4 passes** if no quantity is a changed equilibrium.

**F5, option-off identity:**

- with the tracer off, the code reproduces the previous code bit for bit;
- with it on, the bulk state is bit-identical to the run with it off. The
  tracer is passive.

**F6, water-path topology.** Every transfer has a recorded destination, and
nothing is left to an unnamed residual:

- the sum of the classes equals the layer's vapour at every step, to 1e-12
  relative;
- P_surface,unique ≤ C_precip,unique ≤ C_unique, cumulatively at every step;
- PE = FE·SE to round-off;
- **inventory closure,** each to 1e-12 relative to the throughput:
  - convective-born: Δ(C and R inventory) = C_unique − P_surface,unique
    (all routes);
  - U: Δ(U inventory) = surface evaporation − first condensation
    (convective and large-scale);
- **the recycling-invariance unit test.** A parcel of mass m is taken through
  condense → evaporate → condense → evaporate → condense → form rain →
  evaporate → condense → form rain → reach the surface. It must give
  C_unique = C_precip,unique = P_surface,unique = m and FE = SE = PE = 1,
  however large its gross totals.

## Decision

**If F1-F6 all pass at all three timesteps,** F0 is recorded as "a
numerically acceptable closure candidate awaiting radiative-coupling
validation".

**If any fails,** that gate's failure is recorded with its numbers, and F0 is
not a candidate until a protocol addresses it.

**Reported either way:**

- RH 300-700 and its profile;
- surface precipitation and temperature;
- CAPE, duty cycle, mass transport and subsidence;
- the gross ledger and the unique quantities;
- R_recycle.

## Code, introduced by this protocol only

- **`track_water_provenance`.** It is an option of the water module, which
  owns the class state as a diagnostic partition of the vapour it owns. It
  is also an option of convection, which computes each step's class changes
  and publishes them. The water module applies them with its own processes,
  as it does convection's vapour transport.
- **Accepted configurations:** the carried pathway with
  precipitation_evaporation "rate_law_only" or "none". Others refuse.
- **`mass_flux.provenance_substep`** propagates the classes through one
  substep: the plume, the operator, detrainment, the evaporation of
  detrained cloud, the conversion, the fall and its evaporation, and the
  surface arrival.
- **Tests:**
  - the recycling-invariance parcel;
  - the classes sum to the vapour through a coupled run;
  - the ordering and PE = FE·SE;
  - both inventory closures;
  - bulk identity with the tracer on and off;
  - option-off identity with the previous code.
