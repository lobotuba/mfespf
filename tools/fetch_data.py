"""Fetch the declared external datasets and verify them against their record.

    python tools/fetch_data.py            # fetch what is missing, verify the rest
    python tools/fetch_data.py --verify   # verify only, fetch nothing
    python tools/fetch_data.py --update   # record the checksum of what is on disk

The bytes are not in the repository; ``data/sources.yaml`` is. That separation
is deliberate. The provenance of a dataset -- where it came from, what it is,
what licence it carries, and the checksum that identifies it -- belongs under
version control, because a result is not reproducible without it. The dataset
itself belongs to whoever produced it, and fetching a copy for your own use is
not the same permission as shipping one.

The script refuses to overwrite a file whose checksum differs from the record,
for the same reason ``RunManifest.write`` refuses to replace retained evidence
(A8.4): a dataset that changed under a result is a fact to be noticed, not a
file to be silently replaced. ``--update`` is how you accept a new version, and
it makes the change visible in the diff.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import urllib.error
import urllib.request
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parent.parent
SOURCES = ROOT / "data" / "sources.yaml"
RAW = ROOT / "data" / "raw"
TIMEOUT = 120
USER_AGENT = "mfespf-fetch-data/1.0 (research use; contact the repository owner)"


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_sources() -> dict[str, Any]:
    return yaml.safe_load(SOURCES.read_text(encoding="utf-8"))


def record_checksum(dataset_id: str, checksum: str, retrieved: str) -> None:
    """Update one source's checksum in place, leaving the rest of the file alone.

    A full YAML re-dump would be shorter and would silently discard every
    comment in the file, including the header explaining why the bytes are not
    in the repository and the note listing what could not be fetched. This is a
    provenance record people read, so the prose is part of it, and a tool that
    eats prose to save ten lines of code is a bad trade.
    """
    lines = SOURCES.read_text(encoding="utf-8").splitlines(keepends=True)
    start = next(
        (i for i, line in enumerate(lines) if line.strip() == f"- dataset_id: {dataset_id}"),
        None,
    )
    if start is None:
        raise KeyError(f"no source block for {dataset_id!r} in {SOURCES}")
    end = next(
        (i for i in range(start + 1, len(lines)) if lines[i].lstrip().startswith("- dataset_id:")),
        len(lines),
    )

    replacements = {
        "sha256": f"  sha256: {checksum}\n",
        "retrieved": f"  retrieved: '{retrieved}'\n",
    }
    seen = set()
    for i in range(start, end):
        for key, replacement in replacements.items():
            if lines[i].startswith(f"  {key}:"):
                lines[i] = replacement
                seen.add(key)

    # A key the block never carried is appended inside that block, not at the
    # end of the file, so a hand-added source without a checksum still works.
    insert_at = end
    while insert_at > start and not lines[insert_at - 1].strip():
        insert_at -= 1
    for key in ("sha256", "retrieved"):
        if key not in seen:
            lines.insert(insert_at, replacements[key])
            insert_at += 1

    SOURCES.write_text("".join(lines), encoding="utf-8")


def download(url: str, destination: Path) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
        payload = response.read()
    destination.write_bytes(payload)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="verify only, fetch nothing")
    parser.add_argument(
        "--update",
        action="store_true",
        help="write the checksum of whatever is on disk into data/sources.yaml",
    )
    parser.add_argument("--only", help="operate on one dataset_id")
    args = parser.parse_args(argv)

    sources = load_sources()["sources"]
    if args.only:
        sources = [s for s in sources if s["dataset_id"] == args.only]
        if not sources:
            print(f"no source with dataset_id {args.only!r}", file=sys.stderr)
            return 2

    RAW.mkdir(parents=True, exist_ok=True)
    changed = False
    failures: list[str] = []

    for source in sources:
        dataset_id = source["dataset_id"]
        path = RAW / source["filename"]
        recorded = source.get("sha256") or ""

        if not path.exists():
            if args.verify:
                failures.append(f"{dataset_id}: {path.name} is not present")
                continue
            print(f"fetching {dataset_id} ...", end=" ", flush=True)
            try:
                download(source["url"], path)
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                print("FAILED")
                failures.append(f"{dataset_id}: {exc}")
                continue
            print(f"{path.stat().st_size} bytes")

        actual = sha256_of(path)
        if not recorded or args.update:
            if recorded and recorded != actual:
                print(f"  {dataset_id}: checksum changed {recorded[:16]} -> {actual[:16]}")
            if not args.verify:
                record_checksum(dataset_id, actual, datetime.now(UTC).isoformat(timespec="seconds"))
            changed = True
        elif actual != recorded:
            failures.append(
                f"{dataset_id}: {path.name} does not match its record\n"
                f"    recorded {recorded}\n"
                f"    on disk  {actual}\n"
                f"    The file changed under the record. Re-run with --update to accept the "
                f"new version, which puts the change in the diff where a reader will see it."
            )
        else:
            print(f"{dataset_id}: verified {actual[:16]}")

    if changed and not args.verify:
        print(f"\nchecksums written to {SOURCES.relative_to(ROOT)}")

    if failures:
        print(f"\n{len(failures)} problem(s):", file=sys.stderr)
        for failure in failures:
            print(f"  {failure}", file=sys.stderr)
        return 1

    print(f"\n{len(sources)} dataset(s) present and matching their record.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
