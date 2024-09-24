Release Notes
=============

Below are the release notes for all bandframework releases.

May reference issues on:
https://github.com/bandframework/bandframework/issues

Release 0.4.0
-------------

:Date: MONTH DAY, 2024

New capabilities and notable changes:

- added BAND-compatible PUQ, a parallel package for generating experimental designs tailored for uncertainty quantification at `v0.1.0 <https://github.com/parallelUQ/PUQ/releases/tag/v0.1.0>`_
- added BAND-compatible jitr, a package containing a Lagrange mesh R-matrix solver for parametric reaction model calibration at `v2.0.1 <https://github.com/beykyle/jitr/releases/tag/v2.0.1>`_
- added BAND-compatible nsat, illustrating a Bayesian mixture model approach to quantifying the empirical nuclear saturation point
- updated BAND-compatible parMOO to `v0.4.1 <https://github.com/parmoo/parmoo/releases/tag/v0.4.1>`_, which now includes JIT compilation and automatic differentiation capabilities via `jax`
- updated BAND-compatible rose to `v1.1.3 <https://github.com/bandframework/rose/releases/tag/v1.1.3>`_, includes new features on the backend for performance and greatly expands our test coverage
- updated BAND-compatible SaMBA to `v1.1.0 <https://github.com/asemposki/SAMBA/releases/tag/v1.1.0>`_
- updated BAND-compatible surmise to `v0.3.0 <https://github.com/bandframework/surmise/releases/tag/v0.3.0>`_, which adds coverage and other features and extends testing and documentation

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
