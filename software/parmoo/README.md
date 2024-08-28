# ParMOO in BAND 

ParMOO is a Python library for parallel multiobjective simulation optimization. It seeks to exploit simulation-based structure in objective and constraint functions by building surrogate models of simulation outputs separately from objectives and constraints. 

A [BAND SDK v0.2 Community Policy](/resources/sdkpolicies/bandsdk.md) compatibility documentation for ParMOO is contained in [parmoo-bandsdk.md](/software/parmoo/parmoo-bandsdk.md).

After covering basic installation and usage, here we provide a short tutorial from the [ParMOO Solver Farm](https://github.com/parmoo/parmoo-solver-farm/) that illustrates how ParMOO can be deployed on problems such as the calibration of energy density functionals. 

## ParMOO Installation

Complete installation and testing details for ParMOO are available at the [ParMOO repo](https://github.com/parmoo/parmoo). ParMOO documentation is available on [readthedocs](https://parmoo.readthedocs.io/).

The easiest way to get a minimal installation of ParMOO is via the Python package index, PyPI (commonly called pip):

```
pip install < --user > parmoo
```
where the braces around `< --user >` indicate that the `--user` flag is optional.

For a full installation with all dependencies, including ``libEnsemble``, use:

```
pip install < --user > "parmoo[extras]"
```

You can also clone ParMOO from our GitHub_ and ``pip`` install it
in-place, so that you can easily pull the latest version or checkout
the ``develop`` branch for pre-release features.
On Debian-based systems with a bash shell, this looks like:


```
git clone https://github.com/parmoo/parmoo
cd parmoo
pip install -e .
```

Alternatively, the latest release of ParMOO (including all required and
optional dependencies) can be installed from the ``conda-forge`` channel using:


```
conda install --channel=conda-forge parmoo
```

Before doing so, it is recommended to create a new conda environment using:


```
conda create --name channel-name
conda activate channel-name
```

For additional information, see [ParMOO's installation docs](https://parmoo.readthedocs.io/en/latest/install.html).

## ParMOO References

Please use one or more of the following to cite ParMOO.

The JOSS paper:
```
@article{parmoo,
    author={Chang, Tyler H. and Wild, Stefan M.},
    title={{ParMOO}: A {P}ython Library for Parallel Multiobjective Simulation Optimization},
    journal = {Journal of Open Source Software},
    volume = {8},
    number = {82},
    pages = {4468},
    year = {2023},
    doi = {10.21105/joss.04468}
}
```

ParMOO's online documentation:

```
@techreport{parmoo-docs,
    title       = {{ParMOO}: {P}ython Library for Parallel Multiobjective Simulation Optimization},
    author      = {Chang, Tyler H. and Wild, Stefan M. and Dickinson, Hyrum},
    institution = {Argonne National Laboratory},
    number      = {Version 0.4.0},
    year        = {2024},
    url         = {https://parmoo.readthedocs.io/en/latest}
}
```

Our design principles paper:

```
@techreport{ParMOODesign24,
    title = {Designing a Framework for Solving Multiobjective Simulation Optimization Problems},
    author = {Tyler H. Chang and Stefan M. Wild},
    institution = {arXiv},
    number = {2304.06881},
    year = {2024},
    url = {https://arxiv.org/abs/2304.06881},
}
```


## ParMOO Tutorial: Calibrating a NN Surrogate of the Fayans Energy Density Functional

To try calibrating a Fayans EDF with ParMOO, clone the [Fayans Model Calibration example](https://github.com/parmoo/parmoo-solver-farm/tree/main/fayans-model-calibration-2022) from the [ParMOO Solver Farm](https://github.com/parmoo/parmoo-solver-farm/):
```
git clone https://github.com/parmoo/parmoo-solver-farm
```

From within the subdirectory ``parmoo-solver-farm/fayans-model-calibration-2022``, make sure that all of the dependencies for this example are installed via:
```
python3 -m pip install -r REQUIREMENTS.txt
```

Test installation and execution by running:
```
python3 parmoo_fayans_test.py --comms local --nworkers 4
```

Running the above test produces a CSV file containing ParMOO's final database of design parameter and objective values in `fayans_test_results.csv` and saves a pairwise scatter plot of the Pareto front in `Pareto Front.png`, which should look similar (i.e., small differences will be seen because of the randomization internal to these routines; reproducibility may be obtained on your machine by fixing the numpy random number generator seed) to

![](https://github.com/parmoo/parmoo-solver-farm/blob/main/fayans-model-calibration-2022/Pareto-Front.png)

See the code in [`parmoo_fayans_test.py`](https://github.com/parmoo/parmoo-solver-farm/blob/main/fayans-model-calibration-2022/parmoo_fayans_test.py) for more details on how to reproduce these results, and see the file [`parmoo_fayans_structured_solver.py`](https://github.com/parmoo/parmoo-solver-farm/blob/main/fayans-model-calibration-2022/parmoo_fayans_structured_solver.py) for an example of a production-ready structure-exploiting run.

Details on this problem, which involves a surrogate of the Fayans functional, can be found in Section 5 of the preprint [ArXiv:2304.06881](https://arxiv.org/abs/2304.06881).

## Additional Support and Tutorials

For additional support, see [ParMOO's online docs](https://parmoo.readthedocs.io/en/latest), beginning from the [quickstart guide](https://parmoo.readthedocs.io/en/latest/quickstart.html).
