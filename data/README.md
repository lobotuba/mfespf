# Data

Nothing in this directory is committed except the record of what belongs here.

```bash
python tools/fetch_data.py            # fetch what is missing, verify the rest
python tools/fetch_data.py --verify   # verify only, fetch nothing
python tools/fetch_data.py --update   # accept a new version of a file
```

`sources.yaml` holds the URL, provider, provenance class, licence, citation and
SHA-256 of every declared dataset. `raw/` holds the bytes and is git-ignored.

## Why the bytes are not in the repository

Two things are easy to conflate and are kept apart here.

The **provenance** of a dataset — where it came from, what kind of thing it is,
what licence it carries, and the checksum that identifies it — belongs under
version control, because a result is not reproducible without it. A checksum in
git and a fetch script beside it make a run reproducible from a clean checkout.

The **dataset itself** belongs to whoever produced it. Fetching a copy for your
own research is not the same permission as redistributing one, and the licences
below differ from each other. Keeping `raw/` out of git means the repository
never becomes a distribution channel for someone else's data by accident.

`fetch_data.py` refuses to overwrite a file whose checksum differs from the
record, for the same reason `RunManifest.write` refuses to replace retained
evidence (A8.4): a dataset that changed under a published result is a fact to be
noticed, not a file to be quietly replaced. `--update` is how you accept a new
version, and it puts the change in the diff where a reader will see it.

## What is here

| Dataset | Variable | Coverage | Provenance class |
| --- | --- | --- | --- |
| `HADCRUT5-GLOBAL-MONTHLY` | surface temperature anomaly | 1850–, monthly | observational |
| `BERKELEY-GLOBAL-ANNUAL` | surface temperature anomaly | 1850–, annual | observational |
| `NOAA-OHC-0-700M-ANNUAL` | ocean heat content | 1955–, yearly | observational |
| `NOAA-OHC-0-2000M-PENTADAL` | ocean heat content | 1957–, pentadal | observational |
| `NRLTSI2-ANNUAL` | total solar irradiance | 1610–, annual | proxy |
| `NRLTSI2-MONTHLY` | total solar irradiance | 1882–, monthly | proxy |

`provenance_class` is the framework's classification (INF-002), not the
provider's. NRLTSI2 is **proxy** rather than observational on purpose: before the
satellite era it is a reconstruction driven by sunspot and facular indices, not a
measurement, and a claim resting on the pre-1978 part inherits the
reconstruction's assumptions. The dataset registry refuses an undeclared
provenance class in any evidence role, so this has to be decided rather than
skipped.

## What is missing, and why

**CERES EBAF** is the single most valuable dataset for this framework and is not
here. It reports top-of-atmosphere longwave and shortwave, all-sky and clear-sky
separately, which is exactly how Stage 4 publishes them — a like-for-like
comparison on five fields the model already produces under those names. It needs
a NASA Earthdata login and the CERES ordering tool, so it has to be fetched by
hand. Order EBAF-TOA monthly global means, put the file in `raw/`, add an entry
to `sources.yaml`, and run `--update`.

Absent for the same reason: **CMIP6 solar and volcanic forcing** (input4MIPs,
ESGF authentication) and **HITRAN** line parameters (registration). HITRAN is
what would remove the standing limitation that the radiation module's absorption
coefficients are free parameters rather than measured properties of any gas.

## Two things to decide before using any of this

**Calibration is spent permanently.** The dataset registry enforces INF-003 one
way and at variable granularity: whatever is used to fit a parameter can never
afterwards be independent validation for the same claim. The ocean heat content
series is the live case — use it to fit the two-reservoir exchange coefficient
and it stops being available as evidence about ocean heat uptake. Decide
deliberately rather than finding out when `check_contamination` fires.

**Anomalies are not temperatures.** HadCRUT5 and Berkeley Earth are anomalies
relative to a baseline period; the model produces absolute temperature.
Differencing the model against its own baseline over the same window is a
transformation, and it has to be recorded as one.

## What this framework will show when run against these

The kernel has no CO₂, no aerosols and no volcanic forcing. Driven by solar
irradiance alone it will not reproduce the observed temperature record, and that
is the most useful thing to do with these files rather than a reason to avoid
them. Registered as a claim with its failure criterion declared beforehand,
*"solar variability alone accounts for the observed global-mean temperature
record"* is a falsifiable statement this framework can refute on the record — and
charter §3.5 requires the refutation to stay in the record once it does.
