# Tutorial | 12C(p,gamma)

IMPORTANT: AZURE2 must be available from the command line in the environment
in which BRICK is run.  Currently, AZURE2 is not publicly available. 

There is a link under the login prompt on the [hosting
site](http://azure.nd.edu) where you may request an account. An account will
allow you to download the source code.

Helpful hints:
  (i) The best place to obtain a copy of Minuit2 (necessary for compiling
  AZURE2) seems to be [here](https://github.com/GooFit/Minuit2).
  (ii) GSL is needed for compiling AZURE2 and can be found
  [here](https://www.gnu.org/software/gsl).

For help with the AZURE2 installation, see AZURE2_INSTALL.md.

NOTE: The command-line interface to AZURE2 is currently only available for the
Linux-compatible version. Windows and macOS versions only support the GUI and
are therefore inaccessible to BRICK.

The tutorial is a Jupyter notebook, tutorial/tutorial.ipynb. The goal is to show
the user how to build a Bayesian model to analyze the 12C(p,gamma) capture
reaction data with phenomenological $R$-matrix theory.

