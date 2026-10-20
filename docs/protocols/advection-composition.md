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

## A0 result (2026-10-20)

**The fix: ba42c7c.** The full gate passed: 684 tests, ruff, mypy and
traceability.

- **The single-transport configurations are unchanged.** Over 40 steps, the
  non-tropical two-column and the single tropical column are bit-identical
  to the code before the fix. All 7322 shared outputs match; the only
  differences are the four renamed circulation fields.
- **The old collision is refused,** naming the field and both producers.
- **Both tropical columns now compose with both transports.** Their water
  and convection modules compose with the exchange's rate plus the
  circulation's.

## A1 result (2026-10-20): the rebaseline

**Validity.** It completed 60000 steps.

- **Conservation residual:** 9.1e-16.
- **Holding-relation closure:** 8.2e-15.
- **S1's gates as V0:** G3 transport closure (dse 4.0e-13 W m⁻², vapour
  2.0e-20) and G5 orientation pass. G4 boundedness fails, as it did in V0.
  S1's verdict is "other (reported as found)", as in V0.

**Late means against V0:**

| | V0 (defect) | A1 (fixed) | A1 − V0 |
|---|---:|---:|---:|
| LE rising, W m⁻² | 159.93 | 159.66 | −0.27 |
| LE sinking | 93.46 | 94.71 | +1.25 |
| H rising / sinking | 1.85 / −0.69 | 1.74 / −0.48 | |
| **C_conv rising** | **74.72** | **77.48** | **+2.76** |
| **C_conv sinking** | **19.17** | **14.80** | **−4.37** |
| R_net,sfc rising / sinking | 193.20 / 193.31 | 193.29 / 193.19 | |
| X_rise | 234.65 | 237.14 | +2.49 |
| X_sink | 112.62 | 109.51 | −3.11 |
| **ΔF** | **142.68** | **147.84** | **+5.16** |
| q₁ rising / sinking, kg kg⁻¹ | 0.01630 / 0.01703 | 0.01634 / 0.01699 | |
| q₁,sink − q₁,rise | +0.00073 | +0.00065 | |
| deficit rising / sinking | 0.00639 / 0.00375 | 0.00635 / 0.00380 | |
| lowest-layer RH rising / sinking | 0.726 / 0.816 | 0.727 / 0.815 | |
| free-tropospheric RH (300–700 hPa) rising / sinking | 0.937 / 0.193 | 0.928 / 0.198 | |
| precipitable water rising / sinking, kg m⁻² | 85.16 / 49.19 | 84.94 / 49.64 | |
| sinking ω500, hPa/day | 26.37 | 26.86 | +1.85% |
| sinking mixed-layer top / inversion, hPa | 827.7 / 810.8 | 827.7 / 810.8 | |
| sinking RH jump at the inversion | 0.453 | 0.444 | |
| sinking boundary-layer q / RH | 0.01665 / 0.895 | 0.01664 / 0.891 | |
| band RH 300–700, band precipitable water | 0.606, 68.48 | 0.604, 68.54 | |

**Cross-inversion budget, sinking (V2's table), A1:**

| term | V0 | A1 |
|---|---:|---:|
| O1, kg m⁻² s⁻¹ | 8e-20 | 9e-20 |
| O4 | 7.78e-6 | 7.78e-6 |
| O5 downward | −3.30e-5 | −3.57e-5 |
| O5 net exchange | −2.70e-5 | −2.85e-5 |
| O4 MSE flux, W m⁻² | 19.45 | 19.46 |

**Thresholds regenerated from A1 by the same formulas:**

| quantity | V0 | A1 |
|---|---:|---:|
| G0 | 122.03 | 127.63 W m⁻² |
| substantial exchange | 61.01 | 63.81 W m⁻² |
| ΔF baseline | 142.68 | 147.84 W m⁻² |
| substantial ΔF | 71.34 | 73.92 W m⁻² |
| material transport | 3.737e-6 | 3.787e-6 kg m⁻² s⁻¹ |

**What the fix changed, recorded as found.**

- **The humidity barely moved.** No lowest-layer q or deficit changed by more
  than 1.3 percent, and the reversed near-surface contrast persists
  (+0.00065 against +0.00073).
- **What changed is convection's surface exchange.** Convection now sees the
  circulation's in-step heating: subsidence warming above the sinking
  boundary layer, and ascent cooling in the rising column. The sinking
  column's C_conv fell by 4.37 W m⁻², and the rising column's rose by 2.76.
- **So ΔF widened by 5.16 W m⁻²,** the wrong direction for the S1 problem.
  The fix corrects the coupling; it does not reduce the redistribution.

**A1 replaces V0 as the baseline.** L0 and L0-null are running from A1's
final state.
