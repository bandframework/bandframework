# BAND Framework software
This contains the core software tools for the BAND Framework. The main BAND Framework README is [here](../README.md). 

**Read [CONTRIBUTING](/CONTRIBUTING.rst), and the [bandframework/resources/bandsdk](/resources/sdkpolicies/) 
referenced therein, before attempting to contribute to this directory.**

As of v0.5.0+dev the following packages are present in this directory:
- Bfrescox (v0.0.1-alpha): A Python wrapper for the Frescox coupled reaction channel code.
- jitr (v3.0.1): A Python package containing a Lagrange mesh R-matrix solver for parametric reaction model calibration.
- LCGP (v1.1.1): A Gaussian process surrogate model for emulating stochastic simulation outputs.
- mooGP (v1.0.0): Surrogate model for vector-valued (multi-output) functions with a multi-output orthogonal Gaussian process.
- OpenBT (v1.2.0): C++, Python, and R packages implementing Bayesian tree models, including regression, model mixing, sensitivity analysis and multiobjective optimization.
- parMOO (v0.5.1): A Python library for parallel multiobjective simulation optimization.
- PUQ (v0.1.1): A Python package for generating experimental designs tailored for uncertainty quantification and featuring parallel implementations.
- pybmc (v0.2.4): A Python package for performing Bayesian model combination on various predictive models.
- SmoothEmulator: A simplex sampler, emulator trainer, and MCMC explorer that employs a smooth emulator.
- surmise (v0.4.0): A surrogate model interface for calibration, uncertainty quantification, and sensitivity analysis.
- Taweret (v1.2.0): A Python package containing multiple Bayesian Model Mixing methods.
- rose (v1.1.8): A reduced-order scattering emulator. Note: no longer maintained, see [rose/README.md](https://github.com/bandframework/rose/blob/v1.1.8/README.md) for details.

Applications of these tools to nuclear-physics problems are provided in the ["BAND software uses"](/BANDsoftware_uses) directory. 

Note that the hash for each git submodule contained in this folder is set to the commit of the associated software element of the current BAND 
Framework release. Any serious bug fixes in one of these software elements will occasion a new BAND Framework release. However, some software elements may have version releases that extend functionality in between BAND Framework releases. Please consult the release history of a particular software element to ensure you are using the version appropriate for your research application. 
