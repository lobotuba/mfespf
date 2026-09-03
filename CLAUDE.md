# Working in this repository

MFESPF is a falsifiability framework first and a climate model second. The
specifications in [`specifications/`](specifications/) are normative and govern
every disagreement — charter v0.2 plus Amendment 1, and the v1.0 Interface
Requirements Specification (197 numbered requirements, 189 mandatory).

Read [`docs/architecture.md`](docs/architecture.md) before changing anything in
`coupling/` or `registries/`.

## Non-negotiables

These are enforced in code, and the enforcement is the point of the project. Do
not weaken a check to make a test pass — if a check fires, either the code under
it is wrong or the check's reasoning needs to change, and the second case needs
saying out loud.

- **One owner per prognostic state.** Other modules submit labeled tendencies;
  the owner commits. There is no method on `Module` that writes foreign state.
- **Conservation closes or the run is invalid.** Never widen a tolerance to get
  a green budget. `assert_conservation=False` is for tests that are *about* a
  failure, never for making one go away.
- **Numerical sources stay separate from physical ones.** Anything a limiter,
  filter or clip does to a conserved quantity is a numerical source.
- **No fixed spectra.** Spectral fields travel with their grid. Regridding
  conserves the integral and logs its residual. There is no shortwave/longwave
  split in the radiative transfer -- only radiation per band, with solar and
  thermal sources. Keep it that way; the SW/LW split is a diagnostic formed
  after the solve.
- **Emissivity follows from albedo band by band** (Kirchhoff). Do not add a
  free-standing emissivity parameter to a spectral surface; if a configuration
  wants a reflective-but-black surface, it says so spectrally.
- **A process that moves matter must move its energy in the same interval.**
  The water module owns the latent energy for exactly this reason. Do not add a
  phase change whose energy is applied separately from its mass.
- **A closure that redistributes within a module must sum to zero.** Convective
  adjustment conserves column enthalpy by construction, and the test asserts it
  per step. Do not replace it with a scheme that only conserves on average.
- **Cloud shortwave and longwave effects are published separately, before any
  net.** They are of opposite sign and comparable size; a net alone cannot tell
  a right model from a doubly wrong one (SPC-010, CLD-001).
- **A term that is not physics gets a name, not a hiding place.** A prescribed
  cloud's water enters the budget as `prescribed_cloud_water_source`, the same
  way a limiter's energy enters as a numerical source.
- **Calibration data is never independent validation.** The dataset registry
  enforces it per claim, at variable granularity.
- **Withheld data is unreadable before the claim freezes.** Do not add a
  convenience accessor.
- **A failure criterion is mandatory on every claim.** A claim that cannot fail
  cannot be registered.

## Conventions

- SI internally. Non-SI units convert only at declared boundaries, and the
  conversion is logged. Unknown unit symbols raise — register them in
  `core/units.py` with an explicit dimension and scale.
- Every new equation gets a registry record *before* a module cites it:
  category, mathematical form, assumptions, approximations, validity domain,
  limiting cases, and at least one correctness test.
- Every module contract needs a cost declaration and at least one null,
  limiting-case or disabled-coupling test. The validator refuses it otherwise.
- Override `Module.parameters()` with everything that determines the module's
  numbers. A parameter missing from it makes two different runs hash alike.
- Physics notation wins over PEP 8 naming: `T`, `C`, `S`, `E` are fine (ruff's
  N803/N806 are disabled for this reason).

## Tests

Every test binds to requirement IDs:

```python
@pytest.mark.req("CON-004", "KER-005")
@pytest.mark.stage(1)
def test_matched_transfers_reconcile_to_zero(): ...
```

Citing an ID that does not exist in `requirements/` fails
`tools/traceability.py --check`. Layout:

- `tests/unit/` — one component, no coupled run; analytic identities belong
  here, evaluated directly against the solver rather than through a coupled run
- `tests/integration/` — a coupled configuration end to end
- `tests/acceptance/` — the IRS §22 acceptance criteria, including the failure
  detection ones

Prefer a test that checks physics against an **independently derived** result —
a closed form, a limiting case, an analytic identity — over one that checks the
code against itself.

## Before committing

```bash
python -m ruff check src tests tools && python -m ruff format --check src tests tools
python -m mypy
python -m pytest -q
python tools/traceability.py            # regenerate docs/traceability.md
```

If you amended a specification, regenerate the requirement registry first:
`python tools/extract_requirements.py` (verification method, status and evidence
fields are preserved across regeneration — they are the human record).

## Adding a validation-ladder stage

Stages 5 onward are unimplemented. Follow charter §20 rather than writing the
physics first: define the claim, derive the equation, declare its validity
domain, implement the simplest valid form, verify units and limiting cases,
verify conservation, verify convergence, compare against an independent
process-level observation, couple it, re-run the gates, run a null control,
quantify uncertainty, record the claim status. Preserve negative results —
charter §3.5 requires failed mechanisms and abandoned formulations to stay in
the record when they affect interpretation.
