# Protocol: the cloud-base exchange the observed subcloud layer requires (B0-B2)

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards. The user chose this
step (2026-10-21) on closing `rising-mixing-variants.md`. It is
parameter-free and offline, introduces no new scheme, and changes nothing.

**It follows from:**

- **`rising-mixing-variants.md`:** the subcloud-only limit is not viable.
  Without vapour export through cloud base, the subcloud layer saturates
  within about 30 steps.
- **`rising-mixing-structure.md`:** the rising column's surface-rooted group
  mixes vapour to a uniform q from the surface to 717 hPa on average, at
  every step.

> How large an exchange through cloud base does the observed rising-regime
> subcloud layer require to stay at its observed humidity, and how does the
> exchange the model's full mixing performs compare with it?

## B0: the required exchange (observations and reanalysis)

**The bulk budget of a steady subcloud mixed layer.** Vapour enters at the
surface (E) and by lateral convergence (C), and leaves through cloud base by
an exchange mass flux M between the mixed-layer humidity q_ML and the
humidity just above cloud base, q_+:

M_req = (E + C) / (q_ML − q_+).

This is the standard mixed-layer budget (e.g. Betts 1976, *J. Atmos. Sci.*
33, 1008–1020; Stevens 2006, *Theor. Comput. Fluid Dyn.* 20, 279–304, for
the bulk form). No parameter enters.

**The inputs, all composited on the existing ±20° annual-mean ω500 rising
mask, cos-latitude weighted:**

| input | source | provenance | value already on record |
|---|---|---|---|
| E | OAFlux v3 `lh` 1991–2020, rising regime (`e1.json`), over L = 2.501e6 J kg⁻¹ | observational analysis | 116.97 W m⁻² |
| q_ML | NCEP R1 q at 1000 hPa (`m2_ncep.json`) | reanalysis | 0.017562 |
| q_+ | NCEP R1 q at 925 hPa, the first level above the observed LCL (948 hPa) | reanalysis | 0.013252 |
| C | **bracketed, not observed:** C = 0 for the lower bound; C = the model's own lateral vapour import into its rising subcloud layer (late mean, M0) for the upper bound | — | — |

**So M_req is reported as a range,** [M_req,low, M_req,high].

## B1: the exchange the model performs (offline, M0's window)

**The cloud-base interface:** the model interface nearest NCEP's LCL
(948.0 hPa), from M0's late-mean interface pressures. Frozen as that rule,
not as an index.

**Per step, from M0's saved state and groups,** with D0's validated chain:

- **the gross upward vapour flux by mixing:** water from the layers below
  the interface that the mixing places above it. That is `mix` applied to
  the provisional vapour masked to below the interface, summed above it,
  over Δt;
- **its exchange mass flux:** M_F = that flux / q_below, where q_below is
  the mass-weighted provisional q below the interface;
- **the net export across the interface,** by operator: mixing, the
  circulation's vertical part, and eddy diffusion;
- **the q difference across the interface,** in the committed state;
- **the fraction of the group's mass below the interface exchanged per
  step,** M_F · Δt / m_below. It is reported because full mixing's
  exchange is set by the coupling step, not by a physical rate.

**Also from M0:** the model's lateral vapour import into its rising layers
below the interface, which is C for M_req,high, and E_rise.

## Checks (tolerances scaled to the round-off of what is summed)

| check | rule |
|---|---|
| input | `m0_late.npz` reproduces `r0_late.npz` (recorded, bf9c7d9); the checksums of `e1.json` and `m2_ncep.json` are recorded with the result |
| mixing replay | the masked decomposition sums to the full mixing: mix(below) + mix(above) = mix(all), to ≤ 1e-12 of the column's vapour, every step |
| below-interface budget | over the late window, E + lateral import + net export terms + condensation + storage close for the layers below the interface, to ≤ 1e-12 of the gross terms |

## Classifications, frozen

| class | condition (late means) |
|---|---|
| **Invalid** | a check fails |
| **Full mixing over-exchanges cloud base** | M_F ≥ 10 × M_req,high |
| **Full mixing under-exchanges** | M_F ≤ 0.1 × M_req,low |
| **Comparable** | otherwise |

**Reported always:**

- [M_req,low, M_req,high] and M_F, in kg m⁻² s⁻¹;
- the exchanged fraction per step;
- the model's and NCEP's q differences across cloud base.

## What this cannot separate

- **A bulk mixed-layer budget.** It ignores radiative and precipitation
  effects on the subcloud vapour; there is no condensate in this subcloud
  layer.
- **C is bracketed, not observed.**
- **q_+ is the first reanalysis level above the LCL,** 23 hPa above it,
  not the value exactly at cloud base.
- **M_F is a property of full homogenisation at the coupling step.** It
  would scale with 1/Δt.
- **The references are annual-mean climatologies on a regime mask;** OAFlux
  is ocean-only, while NCEP's mask includes land.

---

## Erratum 1: citations (before B0 reads any data)

**Betts (1976).**

- **As frozen:** *J. Atmos. Sci.* 33, 1008–1020. That page range was
  written from memory and is wrong.
- **As verified:** Betts, A. K., 1976, "Modeling subcloud layer structure
  and interaction with a shallow cumulus layer", *J. Atmos. Sci.* 33,
  **2363–2382**. Betts's own publication list gives the same details.

**Stevens (2006),** "Bulk boundary-layer concepts for simplified models of
tropical dynamics", *Theor. Comput. Fluid Dyn.*

- The paper and its subject are confirmed from its preprint, submitted in
  August 2005, which covers bulk subcloud mixed layers.
- Volume 20, pp. 279–304, is **not verified** against the publisher's
  record.

The method does not depend on either citation: M_req is the subcloud layer's
steady bulk vapour budget, which needs no parameter.

## Results (2026-10-21)

Offline: `scratchpad/b.py` → `b.json`. Input SHA-256:

- `e1.json` 0c38f0ba…5374;
- `m2_ncep.json` 1d7ad0d4…3e3c;
- `m0_late.npz` ee35dbab…2983.

**The cloud-base interface:** index 2, at 945.8 hPa, against NCEP's LCL at
948.0 hPa. Two layers lie below it.

**Checks:**

| check | result | rule |
|---|---|---|
| mixing replay | 3.2e-17 | ≤ 1e-12 |
| below-interface budget | 1.8e-16 | ≤ 1e-12 |

### B0: required exchange

| | value |
|---|---|
| E_obs (OAFlux, rising) | 4.677e-5 kg m⁻² s⁻¹ (116.97 W m⁻²) |
| q_ML (NCEP 1000 hPa) / q_+ (NCEP 925 hPa) | 0.017562 / 0.013252 |
| Δq across cloud base | 0.004310 |
| C, upper bracket (the model's lateral import below the interface) | 4.238e-5 kg m⁻² s⁻¹ |
| **M_req** | **[0.0109, 0.0207] kg m⁻² s⁻¹** |

### B1: the model's exchange (rising column, M0's late window)

| quantity | value |
|---|---|
| gross upward vapour flux by mixing | 2.35e-3 kg m⁻² s⁻¹ |
| **M_F** | **0.147 kg m⁻² s⁻¹** |
| fraction of the below-interface mass exchanged per step | **0.77** |
| q difference across the interface (committed) | 9.7e-6 (NCEP: 0.0043) |
| E_rise | 6.40e-5 kg m⁻² s⁻¹ (1.37 × OAFlux) |

**The below-interface budget, over E_rise:**

| term | / E_rise |
|---|---:|
| export by mixing | 1.049 |
| export by vertical circulation | 0.613 |
| export by diffusion | 0.0002 |
| lateral import | 0.662 |
| condensation | 0 |
| storage | 0.0004 |

**M_F / M_req: 13.6 against the lower bound, 7.1 against the upper.**

### Classification, as frozen

- **Not invalid:** both checks pass.
- **Over-exchanges needs M_F ≥ 10 × M_req,high = 0.207:** M_F = 0.147, so
  no.
- **Under-exchanges:** no.
- **Frozen verdict: Comparable.**

**Recorded as found, beside the class.**

- **The ratio lies at 7–14,** straddling the frozen 10× line depending on
  the lateral bracket.
- **Full mixing exchanges 77% of the subcloud layer's mass with the air
  above, every one-hour step.** That is a property of full homogenisation
  at the coupling step (it scales with 1/Δt), not a physical rate.
- **It holds the model's q difference across cloud base at 1e-5,** against
  the observed 0.0043, a factor of 440.
- **The model's subcloud layer also exports by ascent** (0.61 E_rise),
  which the single bulk M_req includes but M_F does not. The model's total
  export through cloud base is 1.66 E_rise, against E_obs + C.
