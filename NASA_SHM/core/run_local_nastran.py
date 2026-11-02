"""Run a copied case with the local Simcenter Nastran and UCI Enterprise license.

Example (connect UCI VPN first):
    python -m core.run_local_nastran workspace/cases/sample_1

The input model is not changed. The solver's exit code alone is not a success
indicator: Nastran may return zero even after a fatal authorization error.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time
from .paths import require_working_output


def inspect_results(case_dir, job_name="sol144"):
    case_dir = Path(case_dir)
    f06 = case_dir / f"{job_name}.f06"
    log = case_dir / f"{job_name}.log"
    fatal, warnings, subcases = set(), set(), set()
    ended = False
    sections = {"displacements": False, "strains": False, "aerodynamic_monitor_loads": False}
    if f06.exists():
        with f06.open(errors="replace") as stream:
            for line in stream:
                if "FATAL" in line:
                    fatal.add(line.strip())
                if "WARNING MESSAGE" in line:
                    warnings.add(line.strip())
                match = re.search(r"^0\s+SUBCASE\s+(\d+)\s*$", line)
                if match:
                    subcases.add(int(match[1]))
                ended |= "* * * END OF JOB * * *" in line
                sections["displacements"] |= "D I S P L A C E M E N T   V E C T O R" in line
                sections["strains"] |= "S T R A I N S" in line
                sections["aerodynamic_monitor_loads"] |= "A E R O D Y N A M I C   M O N I T O R" in line
    runtime_errors = set()
    elapsed_seconds = None
    analysis_completed = False
    if log.exists():
        with log.open(errors="replace") as stream:
            for line in stream:
                if re.search(r"\b(?:ERROR|FATAL)\b", line):
                    runtime_errors.add(line.strip())
                match = re.search(r"Real:\s+([\d.]+) seconds", line)
                if match:
                    elapsed_seconds = float(match[1])
                analysis_completed |= "Analysis complete  0" in line
    completed = ended and not fatal and analysis_completed
    return {
        "status": ("completed_with_runtime_errors" if runtime_errors else "completed") if completed else "failed_or_incomplete",
        "solver_completed": completed,
        "fatal_messages": sorted(fatal),
        "runtime_errors": sorted(runtime_errors),
        "warning_messages": sorted(warnings),
        "subcases_in_output": sorted(subcases),
        "output_sections": sections,
        "elapsed_seconds": elapsed_seconds,
        "output_files": {p.name: p.stat().st_size for p in case_dir.glob(f"{job_name}.*")
                         if p.suffix.lower() in {".f06", ".f04", ".log", ".op2", ".xdb", ".pch"}},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("--deck", default="sol144.bdf")
    parser.add_argument("--solver", type=Path, default=Path(r"D:\Siemens\Femap 2406\nastran\bin\nastran.exe"))
    parser.add_argument("--license-server", default="29000@license-siemens.eng.uci.edu")
    parser.add_argument("--memory", default="1GB")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--inspect-only", action="store_true")
    args = parser.parse_args()
    case_dir = args.case_dir.resolve(strict=True)
    try:
        require_working_output(case_dir)
    except ValueError as exc:
        parser.error(str(exc))
    deck = (case_dir / args.deck).resolve(strict=True)
    if deck.parent != case_dir:
        parser.error("--deck must name a file directly in case_dir")
    if args.threads < 1:
        parser.error("--threads must be positive")
    invocation = case_dir / "run-invocation.json"
    metadata = json.loads(invocation.read_text(encoding="utf-8")) if args.inspect_only and invocation.exists() else {}
    previous_summary = case_dir / "run-summary.json"
    if args.inspect_only and previous_summary.exists():
        previous = json.loads(previous_summary.read_text(encoding="utf-8"))
        for key in ("launcher_returncode", "fresh_f06"):
            if key in previous:
                metadata[key] = previous[key]
    if not args.inspect_only:
        solver = args.solver.resolve(strict=True)
        scratch = case_dir / "scratch"
        scratch.mkdir(exist_ok=True)
        command = [str(solver), deck.name, f"auth={args.license_server}", "program=NASTRAN",
                   f"memory={args.memory}", f"parallel={args.threads}", f"sdirectory={scratch}",
                   "scr=yes", "notify=no", "old=yes"]
        environment = os.environ.copy()
        environment["SPLM_LICENSE_SERVER"] = args.license_server
        environment["OMP_NUM_THREADS"] = str(args.threads)
        environment["MKL_NUM_THREADS"] = str(args.threads)
        environment.pop("NOLICENSECHECK", None)
        started = time.time()
        metadata = {
            "command": command, "started_unix": started,
            "input_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in case_dir.iterdir() if p.suffix.lower() in {".bdf", ".inc"}},
        }
        (case_dir / "run-invocation.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        with (case_dir / "launcher.stdout.txt").open("wb") as stdout, (case_dir / "launcher.stderr.txt").open("wb") as stderr:
            result = subprocess.run(command, cwd=case_dir, env=environment, stdout=stdout, stderr=stderr,
                                    creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        metadata["launcher_returncode"] = result.returncode
        f06 = deck.with_suffix(".f06")
        metadata["fresh_f06"] = f06.exists() and f06.stat().st_mtime >= started - 2
    summary = inspect_results(case_dir, deck.stem)
    if metadata:
        if metadata.get("fresh_f06") is False:
            summary["status"] = "failed_or_incomplete"
            summary["solver_completed"] = False
        elif metadata.get("launcher_returncode", 0) != 0:
            summary["runtime_errors"].append(f"Launcher returned exit code {metadata['launcher_returncode']}")
            summary["status"] = "completed_with_runtime_errors" if summary["solver_completed"] else "failed_or_incomplete"
    summary.update(metadata)
    (case_dir / "run-summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0 if summary["status"] == "completed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
