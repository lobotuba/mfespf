"""Requirement traceability matrix.

IRS section 24 and ACC-001: every mandatory requirement needs an acceptance
test, inspection, analysis or demonstration record.

Tests bind themselves to requirements with ``@pytest.mark.req("CON-004", ...)``.
This tool reads those markers statically -- it parses the test sources rather
than importing them -- cross-references the registry under ``requirements/``, and
writes ``docs/traceability.md``.

Four failure modes are checked, and they are the ones that matter:

* a test cites a requirement ID that does not exist (a typo silently claiming
  coverage that was never specified);
* a requirement is marked ``verified`` with no evidence behind it;
* a requirement is marked ``carved_out`` with no waiver naming the clause that
  defers it;
* a requirement is marked ``carved_out`` and has a bound test, which is a
  contradiction in one direction or the other.

Coverage of the remaining mandatory requirements is *reported*, not enforced.
Version 1.0 is not finished, and a tool that failed the build for unimplemented
requirements would only encourage marking them done.

Carve-outs (SYS-002)
--------------------

Charter 4.3 defines an architectural carve-out as a reserved and documented
place for a domain, guaranteeing that an implementation may be added and
asserting nothing about whether it is important or active. Charter 18 lists the
fifteen candidate active modules of the first executable kernel and says that
exclusion from it "does not remove the corresponding architectural interface";
IRS 23 repeats it for nine named extensions.

Counting those requirements as unstarted work makes the headline figure mean
less than it should: it reads as a backlog when part of it is a decision. So a
carve-out is a status here, it must cite the clause that defers it, and it is
reported against its own denominator rather than silently inflating the other.

**The rule that decides the hard cases: a carve-out is a domain the kernel does
not use.** A domain it uses without the declared interface is a gap. Gravity
forces the distinction -- the gravity module is as absent from charter 18's list
as the cryosphere is, but the column computes hydrostatic layer mass with it on
every step, out of a bare constant whose own note says the module supersedes it.
GRV-001 is not carved out; it is unimplemented.
"""

from __future__ import annotations

import argparse
import ast
import sys
from collections import defaultdict
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
REQUIREMENTS_DIR = ROOT / "requirements"
TESTS_DIR = ROOT / "tests"
OUTPUT = ROOT / "docs" / "traceability.md"


def load_requirements() -> dict[str, dict]:
    requirements: dict[str, dict] = {}
    for path in sorted(REQUIREMENTS_DIR.glob("req-*.yaml")):
        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        for row in data.get("requirements", []):
            row["group"] = data.get("group", "")
            row["group_title"] = data.get("title", "")
            requirements[row["id"]] = row
    return requirements


def collect_bindings() -> dict[str, list[str]]:
    """Map requirement ID to the tests that claim to verify it."""
    bindings: dict[str, list[str]] = defaultdict(list)
    for path in sorted(TESTS_DIR.rglob("test_*.py")):
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except SyntaxError as exc:  # pragma: no cover - broken test file
            print(f"{path}: {exc}", file=sys.stderr)
            continue
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            for requirement_id in _requirement_marks(node):
                bindings[requirement_id].append(f"{path.relative_to(ROOT).as_posix()}::{node.name}")
    return bindings


def _requirement_marks(node: ast.FunctionDef | ast.AsyncFunctionDef) -> list[str]:
    found: list[str] = []
    for decorator in node.decorator_list:
        if not isinstance(decorator, ast.Call):
            continue
        target = decorator.func
        # Match pytest.mark.req(...) and mark.req(...) alike.
        if not (isinstance(target, ast.Attribute) and target.attr == "req"):
            continue
        for arg in decorator.args:
            if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                found.append(arg.value)
    return found


CARVED_OUT = "carved_out"


def is_carve_out(row: dict) -> bool:
    return row.get("status") == CARVED_OUT


def build_matrix(requirements: dict[str, dict], bindings: dict[str, list[str]]) -> str:
    by_group: dict[str, list[dict]] = defaultdict(list)
    for row in requirements.values():
        by_group[row["group"]].append(row)

    total = len(requirements)
    mandatory = [r for r in requirements.values() if r["strength"] == "mandatory"]
    carved = [r for r in mandatory if is_carve_out(r)]
    expected = [r for r in mandatory if not is_carve_out(r)]
    covered = [r for r in expected if bindings.get(r["id"])]

    lines = [
        "# Requirement traceability matrix",
        "",
        "Generated by `tools/traceability.py`. Do not edit by hand; edit the",
        "registry under `requirements/` or the `@pytest.mark.req` markers in the tests.",
        "",
        f"- requirements: **{total}**",
        f"- mandatory (SHALL / SHALL NOT): **{len(mandatory)}**",
        f"- architectural carve-outs (SYS-002, charter 4.3 and 18): **{len(carved)}**",
        f"- mandatory outside a carve-out: **{len(expected)}**",
        f"- of those, with a bound test: **{len(covered)}** ({len(covered) / len(expected):.0%})"
        if expected
        else "",
        "",
        "Coverage below 100% is expected while the reference kernel is incomplete.",
        "A requirement with no bound test and no waiver is unverified, and this table",
        "is where that shows.",
        "",
        "A carve-out is a reserved place in the architecture for a domain this kernel",
        "does not use, not a task nobody has started; each one names the clause that",
        "defers it. A domain the kernel *does* use without its declared interface is a",
        "gap and is counted as one.",
        "",
    ]

    for group in sorted(by_group):
        rows = sorted(by_group[group], key=lambda r: r["id"])
        title = rows[0].get("group_title") or group
        group_mandatory = [r for r in rows if r["strength"] == "mandatory" and not is_carve_out(r)]
        group_covered = [r for r in group_mandatory if bindings.get(r["id"])]
        share = (
            f"{len(group_covered)}/{len(group_mandatory)}"
            if group_mandatory
            else ("all carved out" if any(is_carve_out(r) for r in rows) else "n/a")
        )
        lines.append(f"## {group} — {title} (mandatory covered: {share})")
        lines.append("")
        lines.append("| ID | Strength | Verification | Status | Evidence |")
        lines.append("| --- | --- | --- | --- | --- |")
        for row in rows:
            tests = bindings.get(row["id"], [])
            evidence = "<br>".join(f"`{t}`" for t in tests) or (
                ", ".join(row.get("evidence") or []) or "—"
            )
            status = row.get("status", "not_started")
            if tests and status == "not_started":
                status = "test_bound"
            lines.append(
                f"| {row['id']} | {row['strength']} | {row.get('verification', 'TEST')} "
                f"| {status} | {evidence} |"
            )
        lines.append("")

    return "\n".join(line for line in lines if line is not None)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="exit non-zero on an integrity error")
    parser.add_argument("--output", default=str(OUTPUT))
    args = parser.parse_args()

    requirements = load_requirements()
    bindings = collect_bindings()

    problems: list[str] = []
    for requirement_id, tests in sorted(bindings.items()):
        if requirement_id not in requirements:
            problems.append(f"unknown requirement {requirement_id!r} cited by {', '.join(tests)}")
    for row in requirements.values():
        if row.get("status") == "verified" and not (
            bindings.get(row["id"]) or row.get("evidence") or row.get("waiver")
        ):
            problems.append(f"{row['id']} is marked verified with no evidence")
        if is_carve_out(row):
            # A carve-out excuses a requirement from the count, so it has to say
            # on whose authority. Without this the status is a way of making the
            # number look better.
            if not row.get("waiver"):
                problems.append(
                    f"{row['id']} is carved out with no waiver naming the clause that defers it"
                )
            if bindings.get(row["id"]):
                problems.append(
                    f"{row['id']} is carved out and has a bound test "
                    f"({', '.join(bindings[row['id']])}); one of the two is wrong"
                )

    matrix = build_matrix(requirements, bindings)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(matrix, encoding="utf-8")

    mandatory = [r for r in requirements.values() if r["strength"] == "mandatory"]
    carved = [r for r in mandatory if is_carve_out(r)]
    expected = [r for r in mandatory if not is_carve_out(r)]
    covered = [r for r in expected if bindings.get(r["id"])]
    print(f"{output.relative_to(ROOT)}: {len(requirements)} requirements, ")
    print(f"  mandatory with a bound test: {len(covered)}/{len(expected)}")
    print(f"  architectural carve-outs, not counted: {len(carved)}")

    if problems:
        print("\ntraceability problems:", file=sys.stderr)
        for problem in problems:
            print(f"  - {problem}", file=sys.stderr)
        if args.check:
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
