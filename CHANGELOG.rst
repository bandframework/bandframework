Release Notes
=============

Below are the release notes for all bandframework releases.

May reference issues on:
https://github.com/bandframework/bandframework/issues

Release 0.5.0
-------------
:Date: December 18, 2025

New organization:

- Refactored the bandframework into BAND software, BAND software uses, and BAND examples.

New capabilities and notable changes:

- added BAND-compatible lcgp at `v0.2.1 <https://github.com/mosesyhc/LCGP/releases/tag/v0.2.1>`_, a tool for latent component Gaussian process emulation for multivariate stochastic simulations.
- added BAND-compatible ModelDiscrepancy at `v1.1.0 <https://github.com/sjaiswal-tifr/ModelDiscrepancy/tree/32d7d2c46009bc8c67a86580a088fb1b747ae029>`_, an example Bayesian framework for model-data comparison that accounts for theoretical uncertainties.
- added BAND-compatible neutron-rich-bmm at `v0.1.0 <https://github.com/asemposki/neutron-rich-bmm/releases/tag/v0.1.0>`_, an example of Gaussian process Bayesian model mixing for the dense matter equation of state.
- added BAND-compatible pybmc at `v0.2.4 <https://github.com/ascsn/pybmc/releases/tag/v0.2.4>`_, a tool for performing Bayesian model combination on various predictive models.
- added BAND-compatible Bfrescox at `v0.0.1-alpha <https://github.com/bandframework/Bfrescox/releases/tag/v0.0.1-alpha>`_, A Python wrapper for the Frescox coupled reaction channels code.
- updated BAND-compatible jitr to `v2.5.1 <https://github.com/beykyle/jitr/releases/tag/v2.5.1>`_, which fixes bugs in calculating some observables; also adds mass tables, a Reaction class, many examples, and other functionality. 
- updated BAND-compatible rose to `v1.1.7 <https://github.com/bandframework/rose/releases/tag/v1.1.7>`_, which includes minor bug fixes.
- updated BAND-compatible SAMBA to `v1.2.0 <https://github.com/asemposki/SAMBA/releases/tag/v1.2.0>`_, which fixes some GP mixing bugs. 
- updated BAND-compatible surmise to `v0.4.0 <https://github.com/bandframework/surmise/releases/tag/v0.4.0>`_, which improves coverage of test suite, integrates Jupyter Book usage examples, and reassigns research (not fully-tested) code.
- updated BAND-compatible Taweret to `v1.2.0 <https://github.com/bandframework/Taweret/releases/tag/v1.2.0>`_, which updates the infrastructure and notebooks, and improves the documentation.
- updated BAND-compatible PUQ to `v0.1.1 <https://github.com/parallelUQ/PUQ/releases/tag/v0.1.1>`_, which extends hetGPy as a base surrogate module and includes two novel sequential design strategies for stochastic simulation models.
- updated BAND example BMEX to `v0.1.4 <https://github.com/massexplorer/bmex-masses/releases/tag/v0.1.4>`_, which includes minor quality of life improvements and bug fixes.
- updated BAND software use case from BAND Camp 2021 to directly include Bfrescox tutorials.

:Known issues:

- We do not yet have a stated policy on the use/documentation of AI tools.
- If one installs all BAND packages in a particular order (and in a single virtual environment, etc.), there may be an incompatibility with the Eigen dependency; fixes are being explored.
- `Notebooks <https://github.com/bandframework/bandframework/tree/develop/BANDsoftware_uses/BFRESCOX>`_  demonstrating the use of the Frescox coupled reaction channels code in BANDsoftware_uses directory need to be updated so that they employ BAND-compatible `Bfrescox <https://github.com/bandframework/Bfrescox/releases/tag/v0.0.1-alpha>`_.


Release 0.4.0
-------------

:Date: October 2, 2024

New capabilities and notable changes:

- added BAND-compatible jitr, a package containing a Lagrange mesh R-matrix solver for parametric reaction model calibration at `v2.0.1 <https://github.com/beykyle/jitr/releases/tag/v2.0.1>`_
- added BAND-compatible `nsat <https://github.com/cdrischler/nuclear_saturation/tree/c4cfa45a1180b2739e217102d7380736d6844a11>`_, illustrating a Bayesian mixture model approach to quantifying the empirical nuclear saturation point
- added BAND-compatible PUQ, a parallel package for generating experimental designs tailored for uncertainty quantification at `v0.1.0 <https://github.com/parallelUQ/PUQ/releases/tag/v0.1.0>`_
- added BAND-compatible `SmoothEmulator </software/SmoothEmulator>`_, a package for building and employing a smooth emulator
- updated BAND-compatible parMOO to `v0.4.1 <https://github.com/parmoo/parmoo/releases/tag/v0.4.1>`_, which now includes JIT compilation and automatic differentiation capabilities via `jax`
- updated BAND-compatible rose to `v1.1.3 <https://github.com/bandframework/rose/releases/tag/v1.1.3>`_, includes new features on the backend for performance and greatly expands our test coverage
- updated BAND-compatible SaMBA to `v1.1.0 <https://github.com/asemposki/SAMBA/releases/tag/v1.1.0>`_
- updated BAND-compatible surmise to `v0.3.0 <https://github.com/bandframework/surmise/releases/tag/v0.3.0>`_, which adds coverage and other features and extends testing and documentation
- updated BAND-compatible Taweret to `v1.1.0 <https://github.com/bandframework/Taweret/releases/tag/v1.1.0>`_, which adds documentation and tutorial notebooks

:Known issues:

- update unit tests for each project under software
- improve documentation
- improve documentation and facilitation of migration from privateband and bandframework branches to bandframework


Release 0.3.0
-------------

:Date: October 10, 2023

New capabilities and notable changes:

- added BAND-compatible BMEX, a web application for exploring nuclear masses and related quantities at `v0.1.1 <https://github.com/massexplorer/bmex-masses/releases/tag/v0.1.1>`_
- added BAND-compatible parMOO, a parallel multiobjective simulation optimization library at `v0.3.1 <https://github.com/parmoo/parmoo/releases/tag/v0.3.1>`_
- added BAND-compatible rose, a reduced-order scattering emulator at `v1.0.0 <https://github.com/bandframework/rose/releases/tag/v1.0.0>`_
- added BAND-compatible Taweret, a package containing multiple Bayesian Model Mixing methods at `v1.0.0 <https://github.com/bandframework/Taweret/releases/tag/v1.0.0>`_
- updated BAND-compatible SaMBA to `v1.0.1 <https://github.com/asemposki/SAMBA/releases/tag/v1.0.1>`_
- updated BAND-compatible surmise to `v0.2.1 <https://github.com/bandframework/surmise/releases/tag/v0.2.1>`_ with new emulation, calibration, and sampling methods
- improved developer guide
- improved overall processes and navigation
- fully resolved continuous integration bug associated with github actions

:Known issues:

- update unit tests for each project under software
- improve documentation
- improve documentation and facilitation of migration from privateband and bandframework branches to bandframework


Release 0.2.0
-------------

:Date: September 23, 2022

Initial release post SDK update.

New capabilities and notable changes:

- updated BAND SDK to v0.2 to reflect state of community testing and documentation
- added BAND-compatible surmise, a surrogate model interface for calibration, uncertainty quantification, and sensitivity analysis
- added BAND-compatible SaMBA, a sandbox for mixing via Bayesian analysis
- added BAND-compatible Bfrescox, a BAND extension of the frescox scattering code for coupled-channels calculations
- added BAND-compatible BRICK, a Bayesian R-matrix inference code kit facilitating extraction of R-matrix parameters from experimental data
- added BAND-compatible QGP_Bayes, a tutorial on the use of JETSCAPE_SIMS tools to infer parameters of the QGP
- added Code of Conduct
- added release process 
- updated overall processes and navigation

:Known issues:

- continuous integration bug associated with github actions
- update unit tests for each project under software
- improve documentation
- add website and other dependencies to release process
- improve documentation and facilitation of migration from privateband and bandframework branches to bandframework

:Desired features:

- enable actions 
