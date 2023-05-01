# ParMOO in BAND 

ParMOO is a Python library for parallel multiobjective simulation optimization. It seeks to exploit simulation-based structure in objective and constraint functions by building surrogate models of simulation outputs separately from objectives and constraints. 

A [BAND SDK v0.2 Community Policy](/resources/sdkpolicies/bandsdk.md) compatibility documentation for ParMOO is contained in [parmoo-bandsdk.md](/software/parmoo/parmoo-bandsdk.md).

After covering basic installation and usage, here we provide a short tutorial from the [ParMOO Solver Farm](https://github.com/parmoo/parmoo-solver-farm/) that illustrates how ParMOO can be deployed on problems such as the calibration of energy density functionals. 

## ParMOO Installation and References

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

ParMOO is also available on conda-forge or it can be installed from source using the `setup.py` file.
For additional information, see [ParMOO's installation docs](https://parmoo.readthedocs.io/en/latest/install.html).

One can cite ParMOO via the paper:
```
@article{parmoo,
    author={Chang, Tyler H. and Wild, Stefan M.},
    title={{ParMOO}: A {P}ython library for parallel multiobjective simulation optimization},
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
    title       = {{ParMOO}: {P}ython library for parallel multiobjective simulation optimization},
    author      = {Chang, Tyler H. and Wild, Stefan M. and Dickinson, Hyrum},
    institution = {Argonne National Laboratory},
    number      = {Version 0.2.2},
    year        = {2023},
    url         = {https://parmoo.readthedocs.io/en/latest}
}
```



## ParMOO Tutorial: Calibrating a NN Surrogate of the Fayans Energy Density Functional

After installing ParMOO, clone the [Fayans Model Calibration example](https://github.com/parmoo/parmoo-solver-farm/tree/main/fayans-model-calibration-2022) from the [ParMOO Solver Farm](https://github.com/parmoo/parmoo-solver-farm/):
```
git clone https://github.com/parmoo/parmoo-solver-farm/fayans-model-calibration-2022
```

From within the cloned directory, install the example via:
```
python3 -m pip install -r REQUIREMENTS.txt
```

Test installation and execution by running:
```
python3 parmoo_fayans_test.py --comms local --nworkers 4
```

Running the above test produces a CSV file containing ParMOO's final database of design parameter and objective values in `fayans_test_results.csv` and saves a pairwise scatter plot of the Pareto front in Pareto Front.png, which should look similar to

![](https://github.com/parmoo/parmoo-solver-farm/blob/main/fayans-model-calibration-2022/Pareto-Front.png)

See the code in `parmoo_fayans_test.py` for more details on how to preproduce these results, and see the file `parmoo_fayans_structured_solver.py` for an example of a production-ready structure-exploiting run.

## Additional Support and Tutorials

For additional support, see [ParMOO's online docs](https://parmoo.readthedocs.io/en/latest), beginning from the [quickstart guide](https://parmoo.readthedocs.io/en/latest/quickstart.html).
