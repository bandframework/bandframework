# bandframework software
This contains the core software tools for the bandframework. 

**Read [bandframework/resources/bandsdk](/resources/sdkpolicies/) and [CONTRIBUTING](CONTRIBUTING.rst) before attempting to contribute to
this directory.**

As of v0.4.0+dev the following packages are present in this directory. 

- surmise ([v0.3.0](https://github.com/bandframework/surmise/releases/tag/v0.3.0 )): A surrogate model interface for calibration, uncertainty quantification, and sensitivity analysis.
- [SmoothEmulator](/software/SmoothEmulator): A simplex sampler, emulator trainer, and MCMC explorer that employs a smooth emulator.
- parMOO ([v0.4.1](https://github.com/parmoo/parmoo/releases/tag/v0.4.1 )): A Python library for parallel multiobjective simulation optimization.
- rose ([v1.1.7](https://github.com/bandframework/rose/releases/tag/v1.1.7 )): A reduced-order scattering emulator.
- Taweret ([v1.1.0](https://github.com/bandframework/Taweret/releases/tag/v1.1.0 )): A Python package containing multiple Bayesian Model Mixing methods.
- PUQ ([v0.1.0](https://github.com/parallelUQ/PUQ/releases/tag/v0.1.0 )): A Python package for generating experimental designs tailored for uncertainty quantification and featuring parallel implementations.
- jitr ([v2.5.1](https://github.com/beykyle/jitr/releases/tag/v2.5.1 )): A Python package containing a Lagrange mesh R-matrix solver for parametric reaction model calibration.

Applications of these tools to nuclear-physics problems are provided in the ["BAND software uses"](/BANDsoftware_uses) directory. 
