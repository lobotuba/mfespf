# Protocol: why dilute CAPE is discontinuous (D)

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards. The user authorised it
(2026-10-13): "1 then 2". It runs after `convective-fine-substep-limit.md`.

**It follows from** `convective-trigger-overshoot.md`:

- **The jumps.** In a post-hoc diagnostic, CAPE(f) along a fractional
  substep fell smoothly, then jumped about 28 J/kg across the trigger within
  1e-9 of f.
- **How often.** 301 of 755 crossings in a 20-day fractional run at 300 s
  were unresolved this way.

**Constraints:**

- Diagnostic only. No model code changes, and no default changes.
- No change to the plume, its CAPE or its termination. Any such change is
  plume physics and the user's decision.

## Read before this was frozen

**What `entraining_plume` computes:**

- **The source:** the layer of largest moist static energy below the source
  floor.
- **The LFC:** the first layer above the source with positive buoyancy.
- **The LNB:** the first layer above the LFC with buoyancy ≤ 0.
- **CAPE_d:** the positive area from the LFC to the LNB, Σ b·dz over the
  layers from the LFC up to the LNB, excluding the LNB itself.
- **The termination ("top"):** the first layer where the integrated buoyancy
  work from the LFC returns to ≤ 0. The plume's mass rises to there, which can
  be above the LNB.

**Erratum to the candidate named in `convective-trigger-overshoot.md`,
both readings kept:**

- **Reading 1 (there):** "the plume terminates at the first layer where its
  integrated buoyancy work turns non-positive", dropping the area above.
- **Reading 2, from the code:**
  - CAPE_d does not depend on the termination.
  - It ends at the first non-positive-buoyancy layer, the LNB.
  - A layer inside the buoyant column whose buoyancy dips to ≤ 0 therefore
    truncates CAPE_d there, and drops all the positive area above it at once.
  - That truncation is a candidate cause. The termination is not.

## Hypotheses

| | cause | what changes across the jump |
|---|---|---|
| **HL, LNB truncation** | An intermediate layer's buoyancy crosses zero, with positive buoyancy above it | The LNB moves down, with the source and LFC unchanged. The CAPE lost equals the positive area between the new and old LNB, to rounding. |
| **HF, LFC change** | The lowest positive layer changes | The LFC moves |
| **HS, source change** | The arg-max of moist static energy switches layer | The source moves |
| **HI, plume inactive** | One side has no LFC | |
| **Other** | none of the above | |

## D1: classify every unresolved crossing

**The run.** The fractional run at 1800/300 s, 20 days, exactly as `of_300`,
in a diagnostic script. The script wraps `fraction_to_trigger` and
`entraining_plume`, recording without changing anything.

**For each unresolved crossing, from the plumes at f_lo and at f_hi = f_lo +
1e-9:**

- the source, the LFC, the LNB (the first layer above the LFC with b ≤ 0,
  read from the plume's buoyancy) and the termination layer;
- CAPE at both;
- the buoyancy of every layer at both;
- the positive area between the two LNBs.

**Each crossing is classified by the first matching row:**

| class | condition |
|---|---|
| HI | either plume is inactive |
| HS | the source differs |
| HF | the LFC differs |
| HL | the LNB differs, and the CAPE difference equals the positive area between the two LNBs within 1e-6 relative |
| Other | anything else |

## D2: where in the column, and how much

**For the HL crossings,** reported:

- the pressure of the truncating layer;
- the buoyancy there on both sides;
- the depth and CAPE of the positive region above it;
- the distribution of jump sizes.

## Verdict, frozen

| verdict | condition |
|---|---|
| **HL supported** | ≥ 90 percent of unresolved crossings are classified HL |
| **Mixed** | otherwise. The class counts are the result. |

**Not concluded either way:**

- whether the CAPE definition or the termination rule is right;
- what an alternative would do.

Those are physics choices for the user. They would be tested in their own
protocol, against the closure's two predictions.

---
