#!/usr/bin/env python
"""Run one CDCC design point, or assemble finished points into the training set.

This is the *expensive* half of the demo and is not executed as part of it: one
design point is roughly 65 core-hours, and the shipped ``data/training.npz``
is the result of 500 of them (~32,500 core-hours on an HPC cluster).  It is
included so the pipeline is complete and auditable -- everything the notebook
consumes can be regenerated from what is in this directory.

    python run_design_point.py --index 42            # one point -> samples/sample_00042.npz
    python run_design_point.py --index 42 --ranks 16 # ... on 16 MPI ranks
    python run_design_point.py --assemble            # samples/ -> ../data/training.npz

Needs ``bfrescoxpro`` with an MPI-enabled frescox (see requirements-model.txt).

Physics
-------
8B + 208Pb -> 7Be + p + 208Pb at 80 MeV/A, continuum-discretized coupled
channels.  Four parameters of the 7Be-p effective interaction are varied; the
central well depth is *not* one of them, because the ground-state ``&Overlap``
carries ``isc=1``, which makes frescox rescale that depth at every sample to
reproduce the 137 keV proton separation energy.

Each design point needs *two* frescox runs, because the two observables want
different continuum discretizations: the angular distribution needs the
continuum out to 10 MeV, the energy distribution a finer grid out to 3 MeV.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"

# ----------------------------------------------------------------- parameters

#: Inference order.  Every theta array in this package uses it.
PARAM_NAMES = ("R_C", "R_WS", "a_WS", "V_so")

#: Interactions held fixed, from the Bfrescox breakup example.
#:
#: Radii here are *reduced*: frescox multiplies every radius of a potential by
#: ``CC = AT**(1/3) + AP**(1/3)`` from that potential's ``type=0`` namelist.
#: For 7Be+208Pb and p+208Pb, AT=208 so CC = 5.925; for p+7Be, AT=1 and AP=0 so
#: CC = 1 and its radii are absolute -- which is why the calibrated R_C and
#: R_WS are absolute fm and the paper's [2, 3] fm ranges make sense directly.

FIXED = {
    "8B_208Pb": {"rC": 2.65},
    "7Be_208Pb": {"rC": 1.3, "V": 114.2, "rV": 1.286, "aV": 0.853,
                  "W": 9.44, "rW": 1.739, "aW": 0.809},
    "p_208Pb": {"rC": 1.3, "V": 34.819, "rV": 1.17, "aV": 0.75,
                "W": 15.34, "rW": 1.32, "aW": 0.601},
}

#: Starting guess for the p+7Be central depth.  frescox refits it (isc=1), so
#: only the search's starting point depends on this.
V_CENTRAL_SEED = 44.675

#: Observable grids.  25 angular points out to 3 degrees and 20 energy bins out
#: to 3 MeV give the paper's 45-dimensional output vector.
N_ANGLES, THETA_MAX_DEG, N_ENERGIES = 25, 3.0, 20


def template_parameters(theta) -> dict[str, float]:
    """Map one parameter vector onto the template's 22 ``@placeholder@`` values.

    R_WS and a_WS each feed *two* placeholders: the paper takes the spin-orbit
    radius and diffuseness to be the same as the central ones.
    """
    theta = np.asarray(theta, dtype=float).ravel()
    if theta.size != 4:
        raise ValueError(f"expected 4 parameters {PARAM_NAMES}, got {theta.size}")
    r_c, r_ws, a_ws, v_so = (float(v) for v in theta)
    by_interaction = {k: dict(v) for k, v in FIXED.items()}
    by_interaction["p_7Be"] = {
        "rC": r_c, "V": V_CENTRAL_SEED, "rV": r_ws, "aV": a_ws,
        "Vso": v_so, "rso": r_ws, "aso": a_ws,
    }
    return {f"{k}_{suffix}": v
            for suffix, params in by_interaction.items()
            for k, v in params.items()}


# ------------------------------------------------------- reading the template

_OVERLAP = re.compile(r"&Overlap\b(.*?)/", re.IGNORECASE | re.DOTALL)


def bin_table(template_path: Path) -> list[dict]:
    """Recover the continuum bin structure by reading the template. In 
    frescox's ``&Overlap`` namelists ``be`` is minus the bin's  midpoint 
    energy and ``er`` is its full width, so a bin spans 
    ``[|be| - |er|/2, |be| + |er|/2]``.  The ground state has no ``er`` and is
    skipped.

    Returns one record per continuum channel, in the order frescox emits them
    into fort.16 -- which is (partial wave, bin) nested, so several records
    share each energy bin.
    """
    text = template_path.read_text()
    rows = []
    for match in _OVERLAP.finditer(text):
        body = " ".join(match.group(1).split())
        fields = dict(re.findall(r"([A-Za-z_][\w()]*)\s*=\s*(\S+)", body))
        if "er" not in fields:
            continue                       # the bound ground state
        centre, width = abs(float(fields["be"])), abs(float(fields["er"]))
        rows.append({"l": int(fields["l"]), "j": float(fields["j"]),
                     "e_lo": centre - width / 2, "e_hi": centre + width / 2})
    if not rows:
        raise ValueError(f"no continuum bins found in {template_path}")
    return rows


# ------------------------------------------------------------- observables

def breakup_channels(results):
    """Every fort.16 channel except elastic, in declaration order."""
    keys = sorted((k for k in results if k != "channel_1"),
                  key=lambda k: int(k.split("_")[1]))
    return [results[k] for k in keys]


def double_differential(results):
    """``(theta_deg, sigma)`` with sigma of shape (n_channels, n_angles), mb/sr."""
    channels = breakup_channels(results)
    theta_deg = np.asarray(channels[0]["Theta_deg"], dtype=float)
    sigma = np.vstack([np.asarray(c["sigma_mb_sr"], dtype=float) for c in channels])
    return theta_deg, sigma


def dsigma_domega(results, angles_deg):
    """Energy-summed breakup angular distribution, **b/sr**.

    Each bin's cross section is already integrated over that bin, so summing
    over channels integrates over the whole continuum.
    """
    theta_deg, sigma = double_differential(results)
    return np.interp(angles_deg, theta_deg, sigma.sum(axis=0)) / 1000.0


def dsigma_de(results, rows):
    """Angle-integrated breakup energy distribution, **mb/MeV**.

    With l <= 3 there are seven (l, j) partial waves per energy bin, so the
    channels must be **summed within each bin** before dividing by the bin
    width.  An s-wave-only calculation has one channel per bin and hides this.
    """
    theta_deg, sigma = double_differential(results)
    if len(rows) != sigma.shape[0]:
        raise ValueError(f"template has {len(rows)} continuum channels but "
                         f"fort.16 has {sigma.shape[0]}")
    theta_rad = np.deg2rad(theta_deg)
    integrated = 2 * np.pi * np.trapezoid(sigma * np.sin(theta_rad), theta_rad, axis=1)

    by_bin: dict[tuple[float, float], float] = {}
    for row, value in zip(rows, integrated):
        key = (row["e_lo"], row["e_hi"])
        by_bin[key] = by_bin.get(key, 0.0) + float(value)

    edges = sorted(by_bin)
    energies = np.array([0.5 * (lo + hi) for lo, hi in edges])
    widths = np.array([hi - lo for lo, hi in edges])
    return energies, np.array([by_bin[k] for k in edges]) / widths


# ------------------------------------------------------------- running frescox

def read_design(path: Path) -> np.ndarray:
    """The (n, 4) design, without needing pandas."""
    rows = np.genfromtxt(path, delimiter=",", names=True)
    return np.column_stack([rows[n] for n in PARAM_NAMES])


def run_one(template: Path, theta, workdir: Path, ranks: int) -> dict:
    """Fill the template, run frescox in ``workdir``, return parsed fort.16."""
    import bfrescoxpro
    from bfrescoxpro import Configuration, parse_fort16

    workdir.mkdir(parents=True, exist_ok=True)
    cfg = Configuration.from_template(template, workdir / "cdcc.in",
                                      template_parameters(theta), overwrite=True)
    # frescox writes fort.* into its cwd, so each run needs its own directory.
    bfrescoxpro.run_simulation(cfg, workdir / "frescox.out", overwrite=True,
                               mpi_setup={"n_processes": ranks}, cwd=workdir)
    fort16 = workdir / "fort.16"
    if not fort16.is_file():
        raise RuntimeError(f"frescox produced no fort.16 in {workdir}")
    return parse_fort16(fort16)


def run_index(index: int, ranks: int, scratch: Path, out_dir: Path,
              keep: bool = False) -> Path:
    design = read_design(DATA / "design.csv")
    if not 0 <= index < len(design):
        raise SystemExit(f"index {index} outside design of {len(design)}")
    theta = design[index]

    # The OpenMP-enabled frescox build refuses to start without this.
    os.environ.setdefault("OMP_NUM_THREADS", "1")

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"sample_{index:05d}.npz"
    if out_path.exists():
        print(f"{out_path} exists; nothing to do")
        return out_path

    print(f"design point {index}: "
          + ", ".join(f"{n}={v:.5f}" for n, v in zip(PARAM_NAMES, theta)))
    # Unique per job so two concurrent runs of the same index cannot collide.
    tag = os.environ.get("SLURM_JOB_ID", str(os.getpid()))
    base = scratch / f"sample_{index:05d}.{tag}"

    angular = run_one(HERE / "cdcc_angular.template", theta, base / "angular", ranks)
    energy = run_one(HERE / "cdcc_energy.template", theta, base / "energy", ranks)

    angles = np.linspace(0.0, THETA_MAX_DEG, N_ANGLES)
    ang = dsigma_domega(angular, angles)
    energies, ene = dsigma_de(energy, bin_table(HERE / "cdcc_energy.template"))
    if ene.size != N_ENERGIES:
        raise RuntimeError(f"expected {N_ENERGIES} energy bins, got {ene.size}")

    np.savez(out_path, theta=theta, angles_deg=angles,
             dsigma_domega_b_per_sr=ang, energies_mev=energies,
             dsigma_de_mb_per_mev=ene, combined=np.concatenate([ang, ene]))
    if not keep:
        shutil.rmtree(base, ignore_errors=True)
    print(f"wrote {out_path}")
    return out_path


def assemble(out_dir: Path, target: Path) -> int:
    """Combine finished samples into the (45, n_theta) matrix the notebook loads."""
    design = read_design(DATA / "design.csv")
    columns, angles, energies, missing = [None] * len(design), None, None, []
    for i in range(len(design)):
        path = out_dir / f"sample_{i:05d}.npz"
        if not path.exists():
            missing.append(i)
            continue
        with np.load(path) as s:
            columns[i] = s["combined"]
            if angles is None:
                angles, energies = s["angles_deg"], s["energies_mev"]
    if angles is None:
        raise SystemExit("no completed samples found")

    d_tot = angles.size + energies.size
    f = np.full((d_tot, len(design)), np.nan)
    for i, col in enumerate(columns):
        if col is not None:
            f[:, i] = col
    # x labels each observable row: (kind, coordinate), kind 0 = angular, 1 = energy.
    x = np.vstack([np.column_stack([np.zeros(angles.size), angles]),
                   np.column_stack([np.ones(energies.size), energies])])
    np.savez(target, theta=design, f=f, x=x,
             angles_deg=angles, energies_mev=energies)
    print(f"wrote {target}: f {f.shape} (observables x samples), "
          f"{len(design) - len(missing)} complete, {len(missing)} missing")
    if missing:
        print(f"  missing indices: {missing[:20]}{'...' if len(missing) > 20 else ''}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--index", type=int, help="design point to run")
    ap.add_argument("--ranks", type=int, default=32, help="MPI ranks for frescox")
    ap.add_argument("--assemble", action="store_true",
                    help="combine finished samples into ../data/training.npz")
    ap.add_argument("--samples", default=str(HERE / "samples"))
    ap.add_argument("--scratch", default=os.environ.get("TMPDIR", "/tmp"))
    ap.add_argument("--out", default=str(DATA / "training.npz"))
    ap.add_argument("--keep-workdir", action="store_true")
    args = ap.parse_args()

    if args.assemble:
        return assemble(Path(args.samples), Path(args.out))
    if args.index is None:
        ap.error("give --index N, or --assemble")
    run_index(args.index, args.ranks, Path(args.scratch), Path(args.samples),
              keep=args.keep_workdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
