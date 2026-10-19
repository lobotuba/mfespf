# Protocol: how much of S1's 143 W m⁻² comes from the shared cloud field? (C0, C1, C2)

**Frozen before any data download, code or run of this protocol.** Results
are appended below a line at the end; nothing above it changes afterwards.
The design is the user's (2026-10-18), recorded here with its implementation,
data source and thresholds fixed. The user approved the source and the
thresholds before freezing.

**It follows from** `coupled-ascent-subsidence.md`. S1, the fixed-contrast
two-column run, was reported as found:

- a dry subsiding branch formed;
- but holding the observed 1.47 K contrast needed ΔF = **+142.7 W m⁻²** of
  surface-energy redistribution (+61.3 into the rising surface, −81.4 out of
  the sinking one);
- that is about 10 times the observed ocean heat-transport contrast (OAFlux
  13.39, NCEP 4.67).

Both columns carried the same prescribed ISCCP tropical-mean cloud field, a
declared confound in exactly that quantity.

**Constraints:**

- diagnostic first, and no new cloud closure;
- no cloud property tuned: amount, optical thickness, top pressure, phase,
  overlap and effective radius stay as the source and the existing module
  give them;
- no default changed;
- F0, interactive-radiation F0 and P1-P3 stay blocked;
- S2 is not opened by this protocol, only made eligible or not (below).

## Question

How much of the 143 W m⁻² implied surface-energy redistribution comes from
applying the same tropical-mean cloud field to both the ascending and the
descending regime? The question is not whether observed clouds can fix the
two-column model. That comes only after this attribution.

## The observational source and the regime composite

**The dataset.** ISCCP H-series Gridded Monthly (HGM).

- **The paper:** Young et al. (2018), *Earth Syst. Sci. Data* 10, 583-593,
  doi:10.5194/essd-10-583-2018.
- **The dataset:** NOAA NCEI, doi:10.7289/V5QZ281S, cited to the version
  used.
- **Years:** 1991 to the end of the H-series record.
- **Daytime cloud-type amounts:** ISCCP types its clouds from visible and
  infrared together, as Rossow and Schiffer's Table 5 does.

**Before use:**

- **size and variables:** the download size and the HGM variable tables (the
  ATBD and NCEI documentation) are checked;
- **checksums:** each file's SHA-256 is recorded in `data/sources.yaml`
  before it is parsed, with `provenance_class: observed`.

**The masks are the existing ones.** The ±20° annual-mean NCEP/NCAR R1
1991-2020 omega500 masks of `ocean-heat-transport-circulation.md`: rising,
omega500 < 0, area fraction 0.558; sinking, the rest. The same cos-latitude
weighting. The ISCCP grid is mapped onto the masks by the same bilinear
interpolation of omega500 the earlier composites used.

**The model's four cloud levels, per regime**
(`cloud_field.ISCCP_TROPICAL_CLOUD_TYPES`):

| level | ISCCP types |
|---|---|
| high | cirrus + cirrostratus |
| deep convective | deep convection |
| middle | altocumulus + altostratus + nimbostratus, liquid and ice |
| low | cumulus + stratocumulus + stratus, liquid and ice |

**For each regime and level:**

- **amount:** the sum of its types' annual-mean daytime amounts;
- **top pressure and visible optical thickness:** the amount-weighted means
  of its types' monthly mean values;
- **ice fraction:** its ice types' amount over its total.

**If HGM does not supply per-type top pressure and optical thickness,** the
case is marked **incomplete**, as frozen, with no value backfilled. A partial
case may be reported only if the missing property can be taken from HGM
itself at a coarser level (for example, per level rather than per type), and
it is labelled as such.

**A check on the composite, reported, not gated.** The same procedure over
the whole ±30° band gives each level's amount, compared with the existing
Table 5 tropical amounts (high 0.211, deep 0.027, middle 0.128, low 0.238;
ISCCP D-series). Differences above 0.05 are flagged, since the H and D series
differ in version and era.

## C0: reproduce S1, saving states

**The run.** S1's configuration exactly (`s1.py`, `coupled-ascent-subsidence.md`):

- surfaces held at 300.26 and 298.79 K;
- area fractions 0.558 and 0.442;
- relaxed spectral WTG with L = 734 km;
- default physics;
- the shared tropical-mean clouds;
- 60000 steps.

**What it saves** for offline radiation: each column's layer temperatures,
vapour and surface temperature every 240 steps over the last fifth, 50
states.

**Identity.** C0 reproduces S1's 60000-step result (`s1_60000.json`) in every
printed value the protocol uses.

## C1: offline cloud-radiative substitution

**On every saved state,** each column's radiation module (the run's own
RRTMG configuration) is solved diagnostically, with nothing fed back. Four
cloud cases:

1. **shared:** the tropical-mean field in both columns. This reproduces C0's
   own fluxes, to 1e-9 W m⁻² relative, the method's check;
2. **regime:** the rising regime's composite in the rising column, and the
   sinking regime's in the sinking column;
3. **clear sky in each column;**
4. *(reported)* each regime's clouds in the other column, to separate cloud
   from atmosphere.

**Per column and case, all-sky and clear-sky:**

- surface net shortwave and net longwave;
- TOA absorbed shortwave and OLR;
- the surface and TOA cloud radiative effects, shortwave and longwave
  **separately**, then net.

**The implied change in the holding transfer.** With every other surface
flux fixed offline, F_hold = −R_net,sfc + LE + H + (convective exchange). So
the change from the shared to the regime case is ΔF_hold = −ΔR_net,sfc in
each column, and:

**ΔF_C1 = ΔF_C0 − (ΔR_net,sfc,rising − ΔR_net,sfc,sinking),** the primary C1
diagnostic, as a mean over the 50 states.

## C2: an interactive fixed-SST run with regime clouds

**Run only if C1 is material** (below). S1's boundaries and coupling,
unchanged. Only the prescribed cloud field changes, to each column's regime
composite.

**The code, introduced only then.** A per-region cloud-field option on
`two_column`, `regional_cloud_types=(rising, sinking)`, passed to each
region's existing `ObservedCloudFieldModule`. None is today's code, bit for
bit, with a unit test of that.

**The run.** As S1: 20000 steps, read over the last fifth, extended in
20000-step increments to at most 60000 under S1's equilibrium criterion.

**Every S1 diagnostic is retained:**

- omega500;
- both columns' RH profiles, band RH and precipitable water;
- the T500 contrast;
- precipitation;
- the WTG heat and moisture exchange;
- F_hold and ΔF;
- atmospheric storage;
- radiation validity;
- boundedness;
- G1-G7.

**Added:** each column's shortwave and longwave cloud effects separately.

## Thresholds, frozen (approved by the user before freezing)

**The reference is ΔF_C0 = +142.7 W m⁻²,** confirmed by C0. ΔF is the
final value: C2's if C2 ran, otherwise C1's.

| term | definition |
|---|---|
| **material** | a reduction ≥ 50 percent, so ΔF ≤ 71.3 W m⁻², with ΔF still > 0 (the ocean delivering into the rising region) |
| **observationally consistent** | 0 < ΔF ≤ 26.8 W m⁻², twice the OAFlux value. The same bound as S1's "energetically unrealistic" verdict. |
| **negligible** | a reduction < 10 percent, so ΔF > 128.4 W m⁻² |
| **wrong sign** | ΔF increases, or becomes ≤ 0 |

## Verdicts, frozen

| verdict | conditions |
|---|---|
| **Cloud confound rejected** | C1: wrong sign or negligible. C2 is not run. |
| **Cloud confound insufficient** | the reduction is ≥ 10 percent, but the final ΔF is > 26.8 W m⁻². That includes C1 between 10 and 50 percent, when C2 is not run. Shared clouds contributed, but did not cause the principal energetic failure. |
| **Atmospheric response offsets the radiative improvement** | C1 is material, but C2's reduction is < 50 percent (ΔF_C2 > 71.3) |
| **Cloud confound supported** | C2 ran and all of the following hold. (1) ΔF_C2 is observationally consistent. (2) Orientation holds (S1's G5). (3) The circulation's free-tropospheric tendency is negative and remains the largest drying term in the descent column. (4) The band is not materially worse than C0: band RH 300-700 up by ≤ 0.03, band precipitable water up by ≤ 5 kg m⁻², and the band's \|N − applied export\| up by ≤ 2 W m⁻². (5) S1's G1, G2, G3 and G4-validity pass. |
| **Anything else** | reported as found |

## S2 eligibility, frozen

**An S2 from the C2 equilibrium becomes eligible only if all of these hold.**
Opening it is then the user's decision.

- **Bounded:** G4 boundedness, with 50-day means drifting by less than 1
  percent.
- **Stored:** the late atmospheric storage rate within ±0.1 W m⁻² in each
  column.
- **Circulating:** descending circulation persists (G5).
- **Drying:** the descent column's drying persists. The circulation is the
  largest free-tropospheric drying term, and descent RH 300-700 is ≤ 0.55 on
  at least 90 percent of late days.
- **Energy:** ΔF is observationally consistent (≤ 26.8 W m⁻²).
- **Valid:** no failure of the radiation-domain validity gate.

**S1's state is not an S2 initializer,** as the user recorded.

## What follows, frozen

| outcome | next |
|---|---|
| clouds explain most of the discrepancy (supported) | re-run the fixed-SST qualification with regime clouds, then consider S2 from that equilibrium |
| insufficient, or offset by the atmosphere | the descending boundary-layer ventilation protocol: why the moist layer is trapped below about 830 hPa, and why the descent column evaporates only 58 percent as much as the rising one |
| rejected | the ventilation protocol, immediately |

**The rising column's saturation stays third.**

## What this cannot separate

- **Observed clouds may not be the right clouds for this model's atmosphere.**
  Regime clouds are prescribed from observation, not predicted. The test
  says how much of ΔF the shared field carries, not whether the model can
  make the right clouds.
- **ISCCP's daytime typing,** its version (H) and its era differ from Table
  5's D-series. The composite check reports that difference.

---

## Result (appended; nothing above this line changed)

### The data and the regime composite

**ISCCP-H HGM basic, v01r00,** 1991-01 to 2017-06.

- **318 monthly files,** 1.3 GB, downloaded from NCEI's `isccp-basic/hgm/`
  directory.
- **Registered before any parse:** every file's SHA-256 is in
  `data/sources.yaml` (`ISCCP-H-HGM-BASIC-YYYY-MM`). The composite script
  re-verifies each checksum before reading.
- **What the basic product supplies:** cloud amount, mean top pressure and
  mean visible optical thickness for all 18 ISCCP types, split by phase. So
  the composite is **complete**, with nothing backfilled.

**Two deviations from the protocol's wording, recorded:**

1. **The registry's provenance class is `observational`.** The protocol wrote
   "observed", which is not a class the registry has.
2. **The 1° rising area fraction is 0.576.** That is on ISCCP's 1° grid,
   against 0.558 on NCEP's 2.5° grid: the same masks (omega500
   interpolated bilinearly), a different grid.

**Each level's amount, top pressure, visible optical thickness and ice
fraction,** annual mean of the monthly climatology, cos-latitude weighted:

| | rising | sinking | ±20° band |
|---|---|---|---|
| high | 0.299, 218 hPa, τ 3.56, ice 0.83 | 0.132, 204 hPa, τ 2.50, ice 0.84 | 0.228, 215 hPa, τ 3.30 |
| deep convective | 0.052, 246 hPa, τ 48.1, ice 0.91 | 0.013, 270 hPa, τ 46.9, ice 0.84 | 0.036, 250 hPa, τ 47.9 |
| middle | 0.151, 559 hPa, τ 6.16, ice 0.00 | 0.107, 584 hPa, τ 5.95, ice 0.00 | 0.132, 568 hPa, τ 6.09 |
| low | 0.152, 809 hPa, τ 4.37, ice 0.00 | **0.267**, 806 hPa, τ 5.06, ice 0.00 | 0.201, 807 hPa, τ 4.76 |

**The shared field S1 used** (ISCCP D-series Table 5, ±30°, with the global
Table 4 top pressures and optical thicknesses), for comparison:

| level | amount | top pressure | τ |
|---|---:|---:|---:|
| high | 0.211 | 267 hPa | 2.2 |
| deep convective | 0.027 | 326 hPa | 35.6 |
| middle | 0.128 | 557 hPa | 4.8 |
| low | 0.238 | 826 hPa | 4.7 |

**The check against Table 5 (reported):** the same procedure over ±30° gives
amounts of 0.207, 0.032, 0.135 and 0.212.

- The differences from Table 5 are −0.004, +0.005, +0.007 and −0.026.
- **None exceeds the 0.05 flag.**

**One artefact, noted.** Zero amounts unpack to −7.4e-7 percent, a
consequence of the float32 scale factor. That makes the low and middle ice
fractions read about −6e-8 instead of 0. It is negligible and left as unpacked.
