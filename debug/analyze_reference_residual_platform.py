"""Validate and plot the reference equation residual experiment."""
from pathlib import Path
import csv
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / "results/reference_residual_platform"


def load(path):
    with path.open() as stream:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]


def main():
    manifest = json.loads((FOLDER / "solver_module_sha256.json").read_text())
    for name, digest in manifest.items():
        assert hashlib.sha256((FOLDER / "source" / name).read_bytes()).hexdigest() == digest
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2), sharey=True)
    report = {"settings": dict(flow="Poiseuille", model="2D3V Shakhov", Pr=2/3,
                              degree=1, triangles=32, velocity=[24,24,24],
                              acceleration=[1,0], stationary_diffuse_walls=True,
                              periodic_streamwise=True, stabilization=[1,1,1,1.875],
                              requested_tolerance=1e-11, maximum_steps=1200),
              "normalization": "Each original/coupled/source norm divided by the same initial zero-state original equation residual (including the force). This is not the stopping criterion.",
              "residual_norm": "sqrt(sum_velocity w * sum_elements R_weak^T M_element^-1 R_weak / total_area)",
              "solver_modules_unchanged": True, "cases": []}
    for ax, delta in zip(axes, [.1, 1, 10]):
        cases = {}
        for mode in ["cis", "gsis"]:
            folder = FOLDER / f"poiseuille_delta{delta:g}_P1" / mode
            assert (folder / "final_status.dat").exists(), folder
            data = load(folder / "equation_residual.csv")
            for row in data:
                assert abs(row["absolute_original"]-row["absolute_source"]) <= row["absolute_coupled"]+1e-12
                assert abs(row["absolute_original"]-row["absolute_coupled"]) <= row["absolute_source"]+1e-12
            if mode == "cis":
                assert max(abs(r["absolute_original"]-r["absolute_coupled"]) for r in data)<1e-12
            last=data[-1]
            tail=[r for r in data if r["step"]>=.75*last["step"]]
            plateau=lambda key: (max(r[key] for r in tail)-min(r[key] for r in tail))/last[key]
            report["cases"].append(dict(delta=delta, mode=mode.upper(), converged=last["macro_update"]<1e-11,
                                        final=last, original_relative_range_last_300_steps=plateau("relative_original"),
                                        coupled_relative_range_last_300_steps=plateau("relative_coupled")))
            cases[mode]=data
        assert abs(cases["cis"][0]["absolute_original"]-cases["gsis"][0]["absolute_original"])<1e-12
        for mode, key, label, style, color in [
            ("cis","relative_original","CIS: original", "-", "#2878b5"),
            ("gsis","relative_original","GSIS: original", "-", "#d95319"),
            ("gsis","relative_coupled","GSIS: coupled", "--", "#7e2f8e"),
            ("gsis","relative_source","GSIS: source gap", ":", "#252525")]:
            data=cases[mode]
            ax.semilogy([r["step"] for r in data],np.maximum([r[key] for r in data],1e-18),
                        style,color=color,label=label,linewidth=1.7)
        ax.set_title(rf"$\delta={delta:g}$")
        ax.set_xlabel("Iteration")
        ax.grid(True,which="both",alpha=.2)
        ax.set_ylim(1e-9,2)
    axes[0].set_ylabel("Equation residual / initial original residual")
    axes[0].legend(fontsize=8,loc="upper right")
    fig.suptitle("Supplied reference solver: Poiseuille, P1, 32 triangles, original 24x24x24 velocity grid",fontsize=11)
    reached=sum(case["converged"] for case in report["cases"])
    fig.text(.5,.015,f"Requested macro-update tolerance: 1e-11; {reached}/6 runs met it. Reaching the step limit is not convergence.",ha="center",fontsize=9)
    fig.tight_layout(rect=(0,.05,1,.94))
    fig.savefig(FOLDER / "residual_platform.png",dpi=180)
    plt.close(fig)
    (FOLDER / "report.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
