# Protocol: a finite, independently constrained cloud-base exchange (M0-M2)

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards. The design is the
user's (2026-10-21), following `cloud-base-export.md` and its disposition
(983dde4). The source search and the user's four decisions are recorded
here, before any code.

**It follows from:**

- **The zero-exchange limit is not viable** (`rising-mixing-variants.md`).
  The subcloud layer saturates, and the cut collapses to layer 0.
- **Full group homogenisation exchanges 77% of the subcloud mass every
  one-hour step.** It holds the cloud-base q jump at about 1e-5, roughly 440
  times below the observed 0.0043, at an exchange mass flux of 0.147 kg m⁻²
  s⁻¹. The independently estimated requirement is 0.011–0.021
  (`cloud-base-export.md`).

**So the missing element is a finite, timestep-independent cloud-base mass
exchange,** neither "mix nothing" nor "homogenise everything".

## The design (the user's)

**A conservative exchange across cloud base of both vapour and the
adjustment's heat variable,** never vapour alone:

- F_q = M_cb (q_sub − q_cloud);
- F_h = M_cb (h_sub − h_cloud).

**The rate comes from an independent observable,** M_cb = ρ_cb a_c w_c:

- a regime-matched observed or independently modelled cloud-base air
  density, active updraft area fraction and in-cloud vertical velocity;
- or a direct observed cloud-base mass flux.

**M_cb is never derived from the model's humidity error,** and never from
M_req (0.011–0.021), which stays an evaluation quantity. It is not
rescaled after M1.

**Stages:**

- **M0:** exact baseline reproduction with the option off.
- **M1:** an offline, no-feedback prediction from saved states.
- **M2:** an interactive fixed-SST run, only if M1 is finite and
  nondegenerate.

**Gates:**

1. conservation;
2. timestep convergence over Δt, Δt/2 and Δt/4;
3. nondegeneracy;
4. a structural comparison, not a fit, against M_req, the observed jump and
   the source's M_cb;
5. the downstream result, with no particular RH required.

**Classes:**

- independently constrained exchange supported;
- transport rate too weak;
- transport rate too strong;
- structurally viable but downstream-insufficient;
- **no usable external constraint:** a regime-matched dataset cannot supply
  M_cb without fitted assumptions. M_req is never substituted.

## The user's decisions, recorded before any code (2026-10-21)

1. **Exchange form: the exact two-reservoir solution.**
   - The surface group splits at cloud base, at
     `rising-mixing-variants.md` erratum 1's fixed point.
   - Each part mixes fully within itself.
   - The two relax toward their joint mean with
     g = 1 − exp(−M_cb Δt (1/m_sub + 1/m_cloud)). Vapour is weighted by
     mass. Heat is in the adjustment's own conserved variable (θ relative to
     the plume target, weighted by C π), so Σ C T stays exact.
   - g = 0 is the split; g = 1 is today's full mixing.
2. **M_cb: both ends of the source's range, each judged separately.** No
   point inside the range is chosen.
3. **Timestep ladder: the current default F is laddered too,** at Δt,
   Δt/2 and Δt/4, so a convergence failure can be attributed.
4. **The rising column: no usable external constraint.** See below.

## The source search (before any code)

**Trade-wind marine, regime-matched to a subsiding trade regime:**

- **Source:** Lamer, K., Kollias, P. & Nuijens, L. (2015), "Observations
  of the variability of shallow trade wind cumulus cloudiness and mass
  flux", *J. Geophys. Res. Atmos.* 120, 6161–6178,
  doi:10.1002/2014JD022950. Two years of ground-based radar and lidar at
  the Barbados Cloud Observatory.
- **What the abstract gives:**
  - nonprecipitating cumuli have a peak cloud fraction of about 5% near a
    mean cloud base of 700 m;
  - the **maximum updraft mass flux is 8–36 g m⁻² s⁻¹ just above cloud
    base;**
  - downdrafts contribute 0–7 g m⁻² s⁻¹.
- **Its averaging basis is inferred, not read.** The full text was not
  reachable (publisher and repository both return 403).
  - As a domain mean with a ≈ 0.05, the in-cloud w would be about
    0.14–0.63 m s⁻¹, which is plausible.
  - As an in-cloud value, it would imply 0.007–0.03 m s⁻¹, which is not.
- **A 2026 preprint at the same site** (egusphere-2026-1974, not yet
  peer-reviewed) finds that updrafts and downdrafts contribute nearly
  equally to the cloud-base mass flux in precipitating clouds.

**Recorded, not used:** this range (0.008–0.036 kg m⁻² s⁻¹) brackets the
budget-derived M_req (0.011–0.021) from `cloud-base-export.md`. The two
share no inputs. Full mixing's 0.147 is 4–18 times the range.

**Ascending regime:**

- **Kumar, Jakob, Protat, Williams & May (2015),** "Mass-flux
  characteristics of tropical cumulus clouds from wind profiler
  observations at Darwin, Australia", *J. Atmos. Sci.* 72, 1837–1855.
- **Not usable here:**
  - it is a land site, under deep convection;
  - its convective area fraction is derived by time–space conversion, not
    measured;
  - no cloud-base mass flux is given in the material reached.

---

## Classification (the user's decision, before M1)

**The rising column: no usable external constraint.**

- The only regime-matched direct cloud-base mass flux found is for trade
  cumulus, a subsiding regime.
- Applying it to the ascending column would be a regime transfer, which
  the user declined.
- No ascending-regime marine source supplying M_cb without fitted
  assumptions was found.
- **M_req is not substituted.**

**M0–M2 are not run.** The rising column's cloud base is where the
donor-layer bias originates (`rising-donor-humidity.md`), so applying the
trade-cumulus value to the sinking column alone would not test the
mechanism this protocol was opened for.

**What would reopen it:** a regime-matched observation of cloud-base mass
flux, or of ρ_cb, a_c and w_c, in the tropical ascending marine regime. The
four decisions above then apply as recorded.

## Final disposition (the user's, 2026-10-21), recorded verbatim

> **No usable external constraint for the rising column.** A cloud-base
> mass-flux range of 0.008–0.036 kg m⁻² s⁻¹ is available for
> nonprecipitating shallow trade cumulus, but not for the ascending tropical
> marine regime responsible for the modeled donor-humidity bias. The
> observed trade-cumulus range overlaps the independently budget-derived
> required exchange of 0.011–0.021 kg m⁻² s⁻¹, while the existing
> full-mixing exchange of 0.147 kg m⁻² s⁻¹ is several times larger. Those
> comparisons are recorded as evidence, not used to set a coefficient. No
> finite-exchange implementation or model run is authorized without a
> regime-matched rising-column constraint.

**The causal chain, as established:**

1. timestep-set full mixing;
2. excessive cloud-base exchange;
3. nearly uniform q through 726–760 hPa;
4. near-saturated donor layers;
5. shallow return flow;
6. 38% of the sinking lowest layer's vapour;
7. the reversed near-surface humidity contrast.

**The correction form, as identified:** a finite two-reservoir conservative
exchange, with g = 1 − exp[−M_cb Δt (1/m_sub + 1/m_cloud)]. What is missing
is not the numerical method but an independently defensible M_cb for the
tropical ascending marine regime.

**The protocol reopens only with a source that gives** M_cb directly, or
enough measured information to calculate M_cb = ρ_cb a_u w_u. It needs:

- cloud-base density ρ_cb;
- the active updraft area fraction a_u;
- the in-cloud cloud-base vertical velocity w_u;
- an averaging definition compatible with a domain-mean model flux;
- a tropical ascending marine regime;
- sufficient distinction between shallow, congestus and deep precipitating
  convection.

**Until then:**

- no midpoint is selected;
- no Barbados value is transferred to the rising column;
- no coefficient is inferred from the desired humidity result;
- no M1 or M2 run proceeds;
- all defaults remain unchanged.

> This closes the investigation at a scientifically defensible boundary: the
> model defect and correction form are identified, but the correction
> amplitude lacks a regime-matched observable.
