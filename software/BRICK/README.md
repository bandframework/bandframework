# Bayesian R-Matrix Inference Code Kit (BRICK)

## About

BRICK is a Python layer around [AZURE2](http://azure.nd.edu). It provides
the functionality to sample R-matrix parameters that are typically fit.

## Workflow

* Set up the analysis in AZURE2. This provides a `.azr` file.
* Instantiate an `AZR` object with the `.azr` file .
* Assign priors.
* Construct a ln(Posterior) function.
* Sample with [emcee](https://emcee.readthedocs.io/en/stable/) (or some other
  sampler).

## Installation

BRICK is hosted on [GitHub](https://github.com/odell/brick). It is also
available on the Python Package Index (PyPI) under the name `brick-james`.

`AZURE2` must be installed for BRICK to work. More details can be found in the
[tutorial README](tutorial/README.md).

## What's Here in BAND?

A tutorial is presented here to lay out the basic functionality of BRICK. It is
a calculation of the 12C(p,gamma) capture reaction.

## Contact Information

BRICK is hosted on [GitHub](https://github.com/odell/brick), so the best way to
contact the development team is to create a new issue or pull request. Let us
know if something isn't working for you or if there's a feature you'd like to
see added.
