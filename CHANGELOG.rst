Release Notes
=============

Below are the release notes for all bandframework releases.

May reference issues on:
https://github.com/bandframework/bandframework/issues

Release 0.3.0
-------------

:Date: October 10, 2023 TBD

Initial release post SDK update.

New capabilities and notable changes:

- added BAND-compatible BMEX, a web application for exploring nuclear masses and related quantities
- added BAND-compatible parMOO, a parallel multiobjective simulation optimization library
- added BAND-compatible rose, a reduced-order scattering emulator
- added BAND-compatible Taweret, a package containing multiple Bayesian Model Mixing methods
- updated BAND-compatible SaMBA
- updated BAND-compatible surmise to `v0.2.0 <https://github.com/bandframework/surmise/releases/tag/v0.2.0>`_ with new emulation, calibration, and sampling methods
- improved developer guide
- improved overall processes and navigation
- github partially fixed continuous integration bug associated with github actions; it is now possible to create and run actions in all bandframework organization repositories except for bandframework/bandframework

:Known issues:

- continuing to experience side effects of the github Actions bug (see `Issue #71 <https://github.com/bandframework/bandframework/issues/71>`_)
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
