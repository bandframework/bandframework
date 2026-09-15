# Emulating and calibrating a nuclear breakup reaction

A self-contained demonstration of Bayesian calibration of an expensive physics simulator,
end to end: **design → simulator outputs → Gaussian-process emulator → validation → MCMC →
posterior predictions → coverage**.

Reproduces Sürer, Nunes, Plumlee & Wild, *Phys. Rev. C* **106**, 024607 (2022), using
[surmise](https://github.com/bandframework/surmise) for the emulation and calibration and
[bfrescoxpro](https://github.com/bandframework/Bfrescox) for the underlying CDCC reaction
model.

## Quick start

Needs Python 3.10–3.13.

**uv** 

```bash
uv sync
OMP_NUM_THREADS=1 uv run jupyter lab breakup_calibration_demo.ipynb
```

Add `--python 3.12` to `uv sync` to force the interpreter the notebook was produced with.

**conda / mamba**

```bash
mamba env create -f environment.yml        # or: conda env create -f environment.yml
conda activate breakup-demo
OMP_NUM_THREADS=1 jupyter lab breakup_calibration_demo.ipynb
```

**pip / venv**

```bash
python3.12 -m venv .venv                   # any of python3.10 .. python3.13
source .venv/bin/activate
pip install -r requirements.txt
OMP_NUM_THREADS=1 jupyter lab breakup_calibration_demo.ipynb
```

On Debian/Ubuntu `python -m venv` fails silently (no `pip`, no `activate`) unless
`python3-venv` is installed: `sudo apt install python3-venv`.

The notebook runs in about five minutes on a laptop — it reads the pre-computed simulator
outputs in `data/` and never invokes the physics code.

`OMP_NUM_THREADS=1` matters. The MCMC does $10^5$ small GP evaluations, and the OpenBLAS that
ships with the numpy wheels (PyPI and conda-forge alike) spreads each one across every core.
On an 8-core machine that made the MCMC cell take 12 minutes instead of 2; the whole notebook
went from 17 minutes to 2.5. If you forget, nothing breaks — it is just slow.

## What is here

```
breakup_calibration_demo.ipynb   the demo
pyproject.toml                   what the notebook needs (source of truth)
uv.lock                          exact versions, for `uv sync`
environment.yml                  the same, for conda / mamba
requirements.txt                 the same, for pip / venv
data/
  training.npz                   500 design points x 45 observables
  mock_data.npz                  simulator at the truth parameters + 10% noise
  design.csv                     the 500-point Latin hypercube design
model/                           the expensive half -- NOT run by the demo
  cdcc_angular.template          frescox input, angular model space (134 states)
  cdcc_energy.template           frescox input, energy model space (141 states)
  run_design_point.py            run one design point, or assemble the .npz
requirements-model.txt           extra dependencies for model/
```

## The problem

$^8$B is a proton halo nucleus — a proton bound to a $^7$Be core by only 137 keV. Fired at
a $^{208}$Pb target at 80 MeV/nucleon it breaks up, and the cross sections depend on four
parameters of the $^7$Be–$p$ interaction: a Coulomb radius, a Woods-Saxon radius and
diffuseness, and a spin-orbit depth. The task is to infer them from the breakup data.

A single CDCC evaluation costs about **65 core-hours**, and MCMC needs $10^5$ of them. So a
500-point design was computed once (~32,500 core-hours on an HPC cluster), a Gaussian
process was fitted to it, and the MCMC runs against the GP.

The "data" are the simulator at known parameter values plus 10% noise, so the posterior can
be checked against a right answer.

## Regenerating the training data

`data/training.npz` is pre-computed. To rebuild it you need a Fortran compiler, an MPI
implementation, and a large machine:

```bash
pip install -r requirements.txt -r requirements-model.txt    # or: uv sync --extra model
BFRESCOX_USE_MPI=enabled BFRESCOX_USE_OPENMP=enabled BFRESCOX_USE_LAPACK=enabled \
  pip install -v "git+https://github.com/bandframework/Bfrescox.git#subdirectory=bfrescoxpro_pypkg"

python model/run_design_point.py --index 0 --ranks 32    # one point, ~2 h on 32 ranks
python model/run_design_point.py --assemble              # samples/ -> data/training.npz
```

In practice this is a SLURM array over `--index 0..499`. Each task is idempotent and uses
its own working directory, so a partial run can be resumed by re-submitting the gaps.
