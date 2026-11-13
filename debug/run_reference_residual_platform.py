"""Run the supplied Fortran solver with read-only equation residual diagnostics.

Solver modules are copied byte-for-byte from Desktop/NonOscillatory. Only the
existing audit driver gains extra monitoring calls. Its original stabilization,
24x24x24 velocity quadrature and Shakhov Pr=2/3 are retained.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import argparse
import csv
import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "reference_runs/nonoscillatory"
OUTPUT = ROOT / "results/reference_residual_platform"
ORIGINAL = Path("/Users/hanchenxi/Desktop/NonOscillatory")


def build():
    source = OUTPUT / "source"
    source.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for path in ORIGINAL.glob("*.f90"):
        if path.name == "Shakhov_2D_3V_DG.f90":
            continue
        shutil.copy2(path, source / path.name)
        manifest[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    for name in ["Pardiso_KLU.c", "Makefile.macos"]:
        shutil.copy2(ORIGINAL / name, source / name)
    # Reuse objects from the already verified reference test build.
    shutil.copytree(REFERENCE / "source/bin_macos", source / "bin_macos", dirs_exist_ok=True)
    driver = (REFERENCE / "source/Shakhov_2D_3V_DG.f90").read_text()
    driver = driver.replace("LOGICAL :: AUDIT_SAVED5", "LOGICAL :: AUDIT_SAVED5\nREAL(KIND=DBL):: AUDIT_PREVIOUS")
    driver = driver.replace("STEP = 1\nDO", """STEP = 0
RESIDUL=1.d0
CALL Audit_Equation_Residual()
AUDIT_PREVIOUS=HUGE(1.d0)
STEP = 1
DO""", 1)
    marker = "\tIF (RESIDUL.LT.TOL) EXIT"
    monitor = """ IF (STEP.LE.2.OR.MOD(STEP,25).EQ.0.OR.RESIDUL.LT.TOL.OR.STEP.GE.TMAX.OR. &
     (RESIDUL.LT.1.d-5.AND.AUDIT_PREVIOUS.GE.1.d-5).OR. &
     (RESIDUL.LT.1.d-7.AND.AUDIT_PREVIOUS.GE.1.d-7).OR. &
     (RESIDUL.LT.1.d-9.AND.AUDIT_PREVIOUS.GE.1.d-9)) CALL Audit_Equation_Residual()
 AUDIT_PREVIOUS=RESIDUL
"""
    assert driver.count(marker) == 1
    driver = driver.replace(marker, monitor + marker, 1)
    include = (ROOT / "debug/reference_equation_residual.inc").read_text()
    driver = driver.replace("END PROGRAM Shakhov_2D_3V_DG", include + "\nEND PROGRAM Shakhov_2D_3V_DG")
    (source / "Shakhov_2D_3V_DG.f90").write_text(driver)
    for name, digest in manifest.items():
        assert hashlib.sha256((source / name).read_bytes()).hexdigest() == digest
    (OUTPUT / "solver_module_sha256.json").write_text(json.dumps(manifest, indent=2) + "\n")
    with (OUTPUT / "build.log").open("w") as log:
        subprocess.run(["make", "-f", "Makefile.macos"], cwd=source, stdout=log, stderr=subprocess.STDOUT, check=True)
    return source / "DGACC_macos"


def run_case(binary, delta, mode, steps):
    spec = importlib.util.spec_from_file_location("reference_channels", REFERENCE / "run_channels.py")
    channels = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(channels)
    folder = channels.prepare_case("poiseuille", delta, 1, 2, 8, mode, 1e-11, OUTPUT)
    control = (folder / "control.in").read_text()
    control, count = re.subn(r"(?m)^(\s*TMAX\s*=)[^\n]*", rf"\g<1> {steps},", control)
    assert count == 1
    (folder / "control.in").write_text(control)
    env = dict(os.environ, REF_MODE=mode, OPENBLAS_NUM_THREADS="1")
    start = time.monotonic()
    with (folder / "run.log").open("w") as log:
        result = subprocess.run([str(binary)], cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT)
    elapsed = time.monotonic() - start
    (folder / "elapsed_seconds.txt").write_text(f"{elapsed:.9g}\n")
    print(f"DONE delta={delta:g} {mode} exit={result.returncode} seconds={elapsed:.1f}", flush=True)
    if result.returncode:
        raise RuntimeError((folder / "run.log").read_text()[-4000:])
    with (folder / "equation_residual.csv").open() as stream:
        rows = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]
    final = rows[-1]
    records = []
    for threshold in [1e-5, 1e-7, 1e-9, 1e-11]:
        hits = [row for row in rows if row["step"] > 0 and row["macro_update"] < threshold]
        if hits:
            records.append(dict(delta=delta, mode=mode, threshold=threshold, reached=True, **hits[0]))
        else:
            records.append(dict(delta=delta, mode=mode, threshold=threshold, reached=False, **final))
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=1200)
    parser.add_argument("--delta", type=float, nargs="+", default=[.1, 1, 10])
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    binary = build()
    print("BUILT reference diagnostic executable", flush=True)
    cases = [(binary, delta, mode, args.steps) for delta in args.delta for mode in ["CIS", "GSIS"]]
    with ThreadPoolExecutor(max_workers=2) as pool:
        batches = list(pool.map(lambda case: run_case(*case), cases))
    rows = [row for batch in batches for row in batch]
    with (OUTPUT / "threshold_summary.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print("SAVED", OUTPUT / "threshold_summary.csv", flush=True)


if __name__ == "__main__":
    main()
