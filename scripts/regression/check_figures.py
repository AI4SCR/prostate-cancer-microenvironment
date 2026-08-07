"""Regression harness: proves a figure/data script's output didn't change.

Usage:
    module load r-light/4.4.1   # only needed for the R entries
    python scripts/regression/check_figures.py                  # fast subset, check against baseline.json
    python scripts/regression/check_figures.py --include-slow   # + the expensive UMAP/sweep scripts
    python scripts/regression/check_figures.py --only figure2_umap.py
    python scripts/regression/check_figures.py --update-baseline --include-slow

Compares rendered output artifacts (PDF/PNG/parquet/CSV), not intercepted
in-memory dataframes -- the terminal representation of "the data that went
into the plot" anyway, and needs no instrumentation of the scripts
themselves. PDFs get their volatile /CreationDate and /ModDate fields
stripped before hashing (matplotlib/ggsurvplot/ComplexHeatmap all embed the
real run timestamp there, which would otherwise make every run "differ"
even with byte-identical content -- confirmed empirically during Phase 1).

Never silences a failure with a per-script tolerance for undiagnosed
non-determinism -- if a script's output isn't reproducible, that's a bug in
the script (see the unseeded-shuffle fix in figure2_umap.py et al.), fixed
at the source, not worked around here.
"""

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

from manifest import MANIFEST, Entry

REPO_ROOT = Path(__file__).resolve().parents[2]
BASELINE_PATH = Path(__file__).resolve().parent / "baseline.json"
_CREATION_DATE_RE = re.compile(rb"/CreationDate\s*\(D:[^)]*\)")
_MOD_DATE_RE = re.compile(rb"/ModDate\s*\(D:[^)]*\)")


def _load_env() -> dict[str, str]:
    """Read DATA_DIR/BASE_DIR/OUTPUT_FIGURES_DIR from .env without a hard
    python-dotenv dependency in this standalone tool (R scripts load .env
    themselves via dotenv::load_dot_env(); Python scripts via
    prostate_cancer.utils; this harness just needs OUTPUT_FIGURES_DIR to
    find what a script wrote)."""
    from dotenv import dotenv_values

    return {**dotenv_values(REPO_ROOT / ".env")}


def output_root_dir(entry: Entry, env: dict[str, str]) -> Path:
    output_figures_dir = Path(env["OUTPUT_FIGURES_DIR"]).expanduser()
    if entry.output_root == "figures":
        return output_figures_dir
    if entry.output_root == "revision":
        return output_figures_dir.parent / "revision"
    raise ValueError(f"unknown output_root {entry.output_root!r}")


def discover_outputs(entry: Entry, env: dict[str, str]) -> list[Path]:
    root = output_root_dir(entry, env)
    paths: set[Path] = set()
    for pattern in entry.outputs:
        paths.update(root.glob(pattern))
    return sorted(p for p in paths if p.is_file())


def hash_file(path: Path) -> str:
    data = path.read_bytes()
    if path.suffix == ".pdf":
        data = _CREATION_DATE_RE.sub(b"/CreationDate (D:STRIPPED)", data)
        data = _MOD_DATE_RE.sub(b"/ModDate (D:STRIPPED)", data)
    return hashlib.sha256(data).hexdigest()


def run_entry(entry: Entry) -> bool:
    """Returns True on success. A script crashing is a failure to report and
    move past, not a reason to abort every other entry in the run."""
    script = REPO_ROOT / entry.script
    if script.suffix == ".py":
        cmd = [sys.executable, str(script)]
    elif script.suffix == ".R":
        rscript = _require("Rscript")
        cmd = [rscript, str(script)]
    else:
        raise ValueError(f"unrecognized script type: {script}")
    print(f"  running {entry.script} ...", flush=True)
    result = subprocess.run(cmd, cwd=REPO_ROOT)
    if result.returncode != 0:
        print(f"  FAIL: {entry.script} exited with code {result.returncode}")
        return False
    return True


def _require(binary: str) -> str:
    import shutil

    path = shutil.which(binary)
    if path is None:
        raise SystemExit(
            f"{binary!r} not found on PATH -- for R entries, run "
            f"`module load r-light/4.4.1` first (see CLAUDE.md)."
        )
    return path


def load_baseline() -> dict[str, dict[str, str]]:
    if not BASELINE_PATH.exists():
        return {}
    return json.loads(BASELINE_PATH.read_text())


def save_baseline(baseline: dict[str, dict[str, str]]) -> None:
    BASELINE_PATH.write_text(json.dumps(baseline, indent=2, sort_keys=True) + "\n")


def select_entries(only: str | None, include_slow: bool) -> list[Entry]:
    if only is not None:
        matches = [e for e in MANIFEST if Path(e.script).name == only or e.script == only]
        if not matches:
            raise SystemExit(f"no manifest entry matches --only {only!r}")
        return matches
    return [e for e in MANIFEST if include_slow or not e.slow]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--include-slow", action="store_true", help="also run the expensive UMAP/sweep scripts")
    parser.add_argument("--only", metavar="SCRIPT", help="run a single entry (filename or repo-relative path)")
    parser.add_argument("--update-baseline", action="store_true", help="overwrite baseline.json with this run's hashes")
    args = parser.parse_args()

    env = _load_env()
    entries = select_entries(args.only, args.include_slow)
    baseline = load_baseline()

    failures: list[str] = []
    for entry in entries:
        print(f"=== {entry.script} ===")
        if not run_entry(entry):
            failures.append(entry.script)
            continue
        outputs = discover_outputs(entry, env)
        if not outputs:
            print(f"  FAIL: no output files matched {entry.outputs} under {output_root_dir(entry, env)}")
            failures.append(entry.script)
            continue

        hashes = {str(p.relative_to(output_root_dir(entry, env))): hash_file(p) for p in outputs}

        if args.update_baseline:
            baseline[entry.script] = hashes
            print(f"  updated baseline ({len(hashes)} files)")
            continue

        expected = baseline.get(entry.script)
        if expected is None:
            print(f"  NEW: not in baseline.json yet ({len(hashes)} files) -- run with --update-baseline to record it")
            continue

        missing = set(expected) - set(hashes)
        extra = set(hashes) - set(expected)
        mismatched = {k for k in expected.keys() & hashes.keys() if expected[k] != hashes[k]}
        if missing or extra or mismatched:
            print(f"  FAIL: {len(missing)} missing, {len(extra)} extra, {len(mismatched)} changed")
            for f in sorted(missing):
                print(f"    missing: {f}")
            for f in sorted(extra):
                print(f"    extra:   {f}")
            for f in sorted(mismatched):
                print(f"    changed: {f}")
            failures.append(entry.script)
        else:
            print(f"  OK ({len(hashes)} files match baseline)")

    if args.update_baseline:
        save_baseline(baseline)
        print(f"\nbaseline.json updated ({len(entries)} entries)")
        return 0

    if failures:
        print(f"\n{len(failures)} script(s) diverged from baseline:")
        for f in failures:
            print(f"  {f}")
        return 1

    print(f"\nall {len(entries)} script(s) match baseline")
    return 0


if __name__ == "__main__":
    sys.exit(main())
