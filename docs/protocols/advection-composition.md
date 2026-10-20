# Protocol: two producers of one advection field (A0-A1)

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards. The user chose
"fix, then rebaseline" (2026-10-20) after L0's smoke test found the defect.
L0 (`boundary-layer-lateral-exchange.md`) runs only after A1, on the
corrected model.

## The defect, as found

**In the tropical two-column configuration, two modules publish the same
fields into each column:**

- the circulation publishes `{region}:vapour_advection_rate` and
  `{region}:advective_heating_rate`;
- each column's `tropical_exchange` (the Brewer–Dobson upwelling and the band
  export) publishes `vapour_advection_rate` and `advective_heating_rate`,
  which the coupler files under the same names.

**The coupler's bus is a dictionary updated in schedule order,** so the later
publisher silently replaces the earlier one. Checked on a 3-step run:

- the water module and the convection module in both columns receive the
  **exchange's** fields, never the circulation's;
- the circulation's transport is still applied, through its own labelled
  tendencies, and every budget closes, which is why nothing flagged it.

**What is lost: the in-step composition.** The water module's saturation
adjustment and the convective adjustment act on a provisional profile that
leaves out this step's circulation transport. They see it a step later, in
the committed state. That breaks the repository's rule that a later module
sees an earlier one's tendency for every quantity.

**Its size, from a 60-step smoke run from V0's final state:**

| missed per step | at most | mean below 600 hPa |
|---|---:|---:|
| circulation heating, rising | 0.157 K | 0.046 K |
| circulation heating, sinking | 0.150 K | 0.064 K |
| circulation vapour, as a fraction of a layer's vapour | 7.4% (rising), 2.5% (sinking) | |

What convection does see, the exchange's heating, is at most 0.016 K per
step.

**Affected:** every tropical two-column run, which covers S1, C0, V0, V1 and
E0–E2. **Unaffected:**

- the non-tropical two-column configuration, where the circulation is the
  only producer;
- the single tropical column, where the exchange is the only producer.

**The tests missed it** because every circulation test uses the non-tropical
two-column configuration.

## A0: the fix (code, before any run)

**1. The coupler refuses two producers of one field.** At configuration
validation, before anything initialises: if two participating modules
declare the same bus name among their produced outputs (after filing under
their regions), the configuration is refused, naming the field and both
producers. Nothing is resolved silently.

**2. Each transport publishes under its own name.**

- The circulation publishes `{region}:circulation_vapour_advection_rate` and
  `{region}:circulation_advective_heating_rate`.
- The exchange keeps `vapour_advection_rate` and `advective_heating_rate`.

**3. The consumers compose with every declared transport.** The water module
and the convection module take a declared tuple of advection sources. They
read each one's rate, declared as a required input, and compose with the sum,
in a fixed order. The sources by configuration:

| configuration | sources |
|---|---|
| single tropical column | the exchange (unchanged) |
| non-tropical two-column | the circulation |
| tropical two-column | the exchange and the circulation |

**4. Nothing else changes.** Each transport's own labelled tendencies, its
transfers and its budgets are unchanged.

**Tests, all before A1:**

- two producers of one field are refused, naming both;
- in the tropical two-column configuration, the water and convection modules
  receive both rates, and compose with their sum;
- **the non-tropical two-column configuration and the single tropical column
  are bit-identical to before the fix,** over a short run, for every shared
  output (one source, so nothing in their arithmetic moves);
- the full gate passes.

## A1: the rebaseline

**The run:** V0's configuration exactly (C0's 60000 steps, the same
instrumentation script), on the fixed code. Read over the last fifth, as V0.

**Validity:** conservation, the holding-relation closure (1e-9), and every S1
gate V0 reported.

**Reported against V0, every quantity V0 reported, and in particular:**

- ΔF, X_rise, X_sink, LE, H and C_conv in both columns;
- the lowest-layer q, RH and surface deficit in both columns;
- the sign of the lowest-layer contrast q₀,sink − q₀,rise;
- the boundary-layer means, the mixed-layer top and the observational
  inversion;
- the cross-inversion budget (V2's table);
- ω500, the free-tropospheric RH and the band precipitable water;
- the thresholds V0 generated, regenerated from A1 by the same formulas.

**No verdict is attached, and nothing is tuned.** A1 replaces V0 as the
baseline for what follows (L0 restarts from A1's final state, with E2-null's
role taken by an A1-configuration restart). The earlier protocols' results
stay as recorded, marked as computed with the defect.

**What A1 can change, recorded before it runs:** only the in-step composition
of convection and condensation with the circulation, not the transport
itself. The direction of any change in the sinking boundary layer is not
predicted.

## The L0 erratum this requires

`boundary-layer-lateral-exchange.md`'s online-split check names "the
published per-layer rate (`{region}:vapour_advection_rate`, and the labelled
circulation tendencies)". In the tropical two-column configuration those two
were not the same field. After A0 the circulation's rate is
`{region}:circulation_vapour_advection_rate`, equal to its labelled tendency.
The erratum is appended to that protocol, recording both readings, before L0
runs.

---
