"""Extract the requirement registry from the Interface Requirements Specification.

The registry under ``requirements/`` is generated from
``specifications/MFESPF_v1.0_Interface_Requirements_Draft.txt`` rather than
retyped, so a requirement cannot exist in the tracker without existing in the
specification, or drift from its wording there.

Run after amending the IRS::

    python tools/extract_requirements.py

Verification method, status and evidence are *not* regenerated: they are the
human judgement the traceability matrix reports (IRS section 24), so existing
values are preserved and only new requirements arrive with placeholders.
"""

from __future__ import annotations

import re
import sys
from collections import defaultdict
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
SPEC = ROOT / "specifications" / "MFESPF_v1.0_Interface_Requirements_Draft.txt"
OUT = ROOT / "requirements"

# A requirement opens a paragraph: three uppercase letters, a hyphen, three
# digits, then its text, which may wrap until a blank line.
_REQ = re.compile(r"^([A-Z]{3})-(\d{3})\s+(.*)$")
_SECTION = re.compile(r"^(\d+)\.\s+(.+)$")

_GROUP_TITLES = {
    "SYS": "System-level requirements",
    "MOD": "Module contract requirements",
    "STA": "State ownership requirements",
    "DAT": "Common field and data requirements",
    "CON": "Conservation and transfer requirements",
    "CAU": "Causality and dependency requirements",
    "HYP": "Hypothesis interface requirements",
    "EQR": "Equation registry requirements",
    "VER": "Interface version requirements",
    "INF": "Information provenance requirements",
    "SPC": "Spectral interface requirements",
    "TIM": "Temporal interface requirements",
    "SPA": "Spatial interface requirements",
    "UNC": "Uncertainty interface requirements",
    "CMP": "Computational resource requirements",
    "AUD": "Run manifest and audit requirements",
    "SOL": "Solar source interface",
    "ORB": "Orbital interface",
    "RAD": "Radiation interface",
    "ATM": "Atmospheric interface",
    "WAT": "Water interface",
    "CLD": "Cloud interface",
    "OCN": "Ocean interface",
    "LND": "Land interface",
    "ICE": "Cryosphere interface",
    "CHM": "Chemistry interface",
    "BIO": "Biosphere and biogeochemistry interface",
    "GRV": "Gravity interface",
    "PAR": "Particle and quantum interface",
    "HUM": "Human activity interface",
    "CLM": "Claim requirements",
    "EXP": "Experiment requirements",
    "VAL": "Validation requirements",
    "KER": "Version 1.0 reference-kernel requirements",
    "ACC": "Version 1.0 acceptance criteria",
}


def classify(text: str) -> str:
    """Normative strength, from the language the requirement actually uses."""
    upper = text.upper()
    if "SHALL NOT" in upper or "SHALL" in upper:
        return "mandatory"
    if "SHOULD" in upper:
        return "recommended"
    if "MAY" in upper:
        return "optional"
    return "unclassified"


def parse(spec_text: str) -> dict[str, list[dict[str, str]]]:
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    lines = spec_text.splitlines()
    section = ""
    i = 0
    while i < len(lines):
        line = lines[i]
        section_match = _SECTION.match(line.strip())
        if section_match and line.strip() == line and not _REQ.match(line):
            section = section_match.group(2).strip()
        match = _REQ.match(line)
        if match:
            prefix, number, first = match.groups()
            body = [first.strip()]
            i += 1
            while i < len(lines) and lines[i].strip() and not _REQ.match(lines[i]):
                if lines[i].startswith("==="):
                    break
                body.append(lines[i].strip())
                i += 1
            text = " ".join(part for part in body if part)
            groups[prefix].append(
                {
                    "id": f"{prefix}-{number}",
                    "text": re.sub(r"\s+", " ", text).strip(),
                    "strength": classify(text),
                    "section": section,
                }
            )
            continue
        i += 1
    return groups


def load_existing(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    return {r["id"]: r for r in data.get("requirements", [])}


def main() -> int:
    if not SPEC.exists():
        print(f"specification not found: {SPEC}", file=sys.stderr)
        return 1
    groups = parse(SPEC.read_text(encoding="utf-8"))
    OUT.mkdir(exist_ok=True)

    total = 0
    for prefix, requirements in sorted(groups.items()):
        # The "req-" prefix is not decoration: an unprefixed "con.yaml" is the
        # reserved Windows device name CON, and git refuses to index it.
        path = OUT / f"req-{prefix.lower()}.yaml"
        existing = load_existing(path)
        rows = []
        for req in requirements:
            prior = existing.get(req["id"], {})
            rows.append(
                {
                    "id": req["id"],
                    "text": req["text"],
                    "strength": req["strength"],
                    "section": req["section"],
                    # Preserved across regeneration: these are the human record.
                    "verification": prior.get("verification", "TEST"),
                    "status": prior.get("status", "not_started"),
                    "design_component": prior.get("design_component", ""),
                    "implementation": prior.get("implementation", ""),
                    "evidence": prior.get("evidence", []),
                    "waiver": prior.get("waiver", ""),
                }
            )
        payload = {
            "group": prefix,
            "title": _GROUP_TITLES.get(prefix, prefix),
            "source": "MFESPF-IRS-001 v1.0 Draft",
            "generated_by": "tools/extract_requirements.py",
            "requirements": rows,
        }
        path.write_text(
            yaml.safe_dump(payload, sort_keys=False, allow_unicode=True, width=100),
            encoding="utf-8",
        )
        total += len(rows)
        print(f"{path.relative_to(ROOT)}: {len(rows)} requirements")

    print(f"\n{total} requirements across {len(groups)} groups")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
