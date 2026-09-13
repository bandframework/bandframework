# pybmc in BAND 

pybmc is a Python package for performing Bayesian Model Combination (BMC) on various predictive models. It provides tools for data handling, orthogonalization, Gibbs sampling, and prediction with uncertainty quantification.

A [BAND SDK v0.2 Community Policy](/resources/sdkpolicies/bandsdk.md) compatibility documentation for pybmc is contained in [pybmc-bandsdk.md](/software/pybmc/pybmc-bandsdk.md).

The version of pybmc adopted for the BAND Framework is v0.4.1, found [here](https://github.com/ascsn/pybmc/tree/v0.3.0).

## pybmc Installation

Complete installation and testing details for pybmc are available at the [pybmc repo](https://github.com/ascsn/pybmc). pybmc documentation is available [here](https://ascsn.github.io/pybmc/).

The easiest way to get a minimal installation of pybmc is via the Python package index, PyPI (commonly called pip):

```
pip install < --user > pybmc
```
where the braces around `< --user >` indicate that the `--user` flag is optional.


You can also clone pybmc from our GitHub and ``pip`` install it
in-place, so that you can easily pull the latest version or checkout
the ``develop`` branch for pre-release features.
On Debian-based systems with a bash shell, this looks like:

```
git clone https://github.com/ascsn/pybmc
cd pybmc
pip install -e .
```

## pybmc Examples

For a simple example of the package in action, check the Usage section of the documentation [here](https://ascsn.github.io/pybmc/usage/#1-load-and-prepare-data). For an interactive notebook of the same procedure, you can load [this file](https://github.com/ascsn/pybmc/blob/v0.3.0/pybmc/test.ipynb) locally or in any Jupyter environment.
