# Protocol: RRTMG's domain, its negative optical depths, and how a refusal is reported

**Frozen before the code.** Results are appended below a line at the end;
nothing above it changes afterwards. This is bookkeeping: what is evaluated,
and how a failure is named. It adds no physics.

## Why

`ocean-heat-transport-circulation.md` N0 and W1 ran away. Each ended 45000
or more steps in, with the solver refusing "cloud optical depths cannot be
negative", recorded as a failed conservation gate. Each part of that report
was wrong or late:

1. **Nothing was a cloud.** The array the solver checks holds every
   absorber, gases included. The negative was RRTMG's.
2. **Nothing failed to conserve.** This is the open `docs/todo.md` item: a
   mid-run refusal is filed as `conservation_reconciliation`.
3. **The column had left RRTMG's tables long before.**
   `EQ-RRTMG-LW-GAS-OPTICS` declares T 160-340 K. AER's k-tables are
   tabulated at the reference profile T_ref(p) and ±15 and ±30 K, and `setcoef`
   extrapolates linearly outside that. N0's tail lay up to 50.5 K outside it
   (`ocean-heat-transport-circulation.md`). The gate could not see this.

**Measured before freezing:** the port yields small negative optical depths
even on RRTMG's own reference profile.

- **An example:** band 7, 880 hPa, 288.6 K, RH 0.5, τ = -0.003.
- **The cause:** the quadratic interpolation in the binary species parameter
  for water-rich mixtures, which is AER's design.
- **AER's handling:** its transfer, `rtrnmr`, sets any negative layer optical
  depth to zero (`if (odepth .lt. 0.0) odepth = 0.0`).
- **This port:** the offline `rtrnmr` already does the same, which is part of
  why P0 matches AER to 5e-5. This model's solver refuses instead.
- **The coupled runs:** R3 and R4 never produced one, or the solver would
  have refused it.

## The three changes

1. **RRTMG's zero-floor, ported.**
   - `rrtmg_optical_depth` sets negative g-point optical depths to zero,
     exactly where `rtrnmr` does.
   - It is not hidden. The radiation module publishes, each step, the count
     of layer-g-point pairs floored and the most negative value before
     flooring.
   - The equation record states it as part of the scheme.
   - It is not a clip of a conserved quantity. An optical depth is not a
     state, and the floor is the published scheme's own.
2. **RRTMG's domain is its tables.**
   - The radiation module publishes, when RRTMG is on, each layer's departure
     from RRTMG's reference temperature at its pressure, T - T_ref(p), with
     T_ref interpolated in log pressure as `setcoef` does.
   - `EQ-RRTMG-LW-GAS-OPTICS` gains the range -30 to +30 K on that departure,
     bound to the field, so the validity gate judges every step.
   - The T 160-340 K range stays: it is the Planck tables' range.
3. **A refusal is named.**
   - When a module raises during a run, the driver records the gate
     `run_completed` as failed, with the error, the module and the step.
   - Only a `ConservationError` fails `conservation_reconciliation`.
   - The solver's message says "absorption optical depths", not "cloud".

## The tests

- **T1:** on RRTMG's reference profile at RH 0.5, `rrtmg_optical_depth`
  returns no negative. Its published floor count is positive and equals the
  number of negatives `rrtmg_gas_optics` gives on the same state.
- **T2:** the published departure equals T - T_ref(p) on a hand-built state.
  A column made 35 K warmer than T_ref in one layer fails the validity gate,
  by name. One within 30 K passes.
- **T3:** a module raising `ConfigurationError` on step 3 leaves a manifest
  whose failed gate is `run_completed`, naming the module, with
  `conservation_reconciliation` not failed. A `ConservationError` still fails
  `conservation_reconciliation`.
- **T4: nothing recorded moves.**
  - P0 and P2 pass unchanged.
  - A 200-step R4-configuration run gives fluxes identical to the code before
    this change, to the bit, with a floor count of zero.

**The full gate passes.**

**Reported, not judged:**

- R3's recorded profile exceeds +30 K in one layer, by 0.3 K at 1.1 hPa. R3
  is not rerun: it was not adopted, and the record stands as it was made.
- R4's largest departure is 27.6 K.
- N0's tail is 50.5 K.

---

## Result (appended; nothing above this line changed)

**The three changes are built, and every test passes.**

- **T1:** on RRTMG's reference profile at RH 0.5, below 100 hPa, nothing
  remains negative after the floor. The floor count equals the raw optics'
  negatives on the grid's g-points, and is positive. The most negative value
  matches to 1e-12.
- **T2:**
  - the published departure equals T - T_ref(p) to 1e-12;
  - a layer started 35 K above the reference fails the validity gate, naming
    `dT_ref` and RRTMG;
  - 25 K above passes.
- **T3:**
  - a `ConfigurationError` on the radiation module's third step fails
    `run_completed` with "ConfigurationError raised by radiation at step 2",
    and `conservation_reconciliation` is not failed;
  - a `ConservationError` still fails `conservation_reconciliation`.
- **T4:**
  - 200 steps of the R4 configuration give all 12400 field values identical
    to the code before the change, to the bit, with masked NaNs compared as
    equal;
  - the floor count is 0 over all 200 steps, and the largest departure is
    27.5 K;
  - P0 and P2 pass unchanged.
- **The full gate passes.**

**Reported.** R3's one layer at +30.3 K (1.1 hPa) and N0's 50.5 K stand as
recorded. Neither is rerun.

**`docs/todo.md`.** The item "A configuration refusal is reported as a
failed conservation gate" is resolved and removed.
