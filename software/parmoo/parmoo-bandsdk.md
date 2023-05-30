# BAND SDK v0.2 Community Policy Compatibility for ParMOO


> This document summarizes the efforts of current and future BAND member packages to achieve compatibility with the BAND SDK community policies.  Additional details on the BAND SDK are available [here](/resources/sdkpolicies/bandsdk.md) and should be considered when filling out this form. The most recent copy of this template exists [here](/resources/sdkpolicies/template.md).
>
> This file should filled out and placed in the directory in the `bandframework` repository representing the software name appended by `bandsdk`.  For example, if you have a software `foo`, the compatibility file should be named `foobandsdk.md` and placed in the directory housing the software in the `bandframework` repository. No open source code can be included without this file.
>
> All code included in this repository will be open source.  If a piece of code does not contain a open-source LICENSE file as mentioned in the requirements below, then it will be automatically licensed as described in the LICENSE file in the root directory of the bandframework repository.
>
> Please provide information on your compatibility status for each mandatory policy and, if possible, also for recommended policies. If you are not compatible, state what is lacking and what are your plans on how to achieve compliance. For current BAND SDK packages: If you were not fully compatible at some point, please describe the steps you undertook to fulfill the policy. This information will be helpful for future BAND member packages.
>
> To suggest changes to these requirements or obtain more information, please contact [BAND](https://bandframework.github.io/team).


**Website:** https://parmoo.readthedocs.io \
**Contact:** parmoo@mcs.anl.gov \
**Icon:** https://raw.githubusercontent.com/parmoo/parmoo/c4de84f3564245d12624081570b6c58923e56dec/docs/img/logo-ParMOO.svg \
**Description:** ParMOO is a Python package for parallel multiobjective optimization that seeks to exploit simulation-based structure in objective and constraint functions.


### Mandatory Policies

**BAND SDK**
| # | Policy                 |Support| Notes                   |
|---|-----------------------|-------|-------------------------|
| 1. | Support BAND community GNU Autoconf, CMake, or other build options. |Full| ParMOO is a Python package and provides a setup.py file for installation. This is compatible with Python's built-in installation feature (``python setup.py install``) and with the pip installer. ParMOO is also available via [Conda-forge](https://anaconda.org/conda-forge/parmoo). GNU Autoconf or CMake are unsuitable for a Python package. |
| 2. | Have a README file in the top directory that states a specific set of testing procedures for a user to verify the software was installed and run correctly. | Full| In addition to a README, ParMOO's installation can be tested via pytest by running ``python3 setup.py test``.|
| 3. | Provide a documented, reliable way to contact the development team. |Full| The ParMOO team can be contacted through the public [issues page on GitHub](https://github.com/parmoo/parmoo/issues) or via an e-mail to [parmoo@mcs.anl.gov](parmoo@mcs.anl.gov).|
| 4. | Come with an open-source license |Full| ParMOO uses the BSD 3-Clause license. |
| 5. | Provide a runtime API to return the current version number of the software. |Full| The version can be returned within Python via: `parmoo.__version__`.|
| 6. | Provide a BAND team-accessible repository. |Full| https://github.com/parmoo/parmoo |
| 7. | Must allow installing, building, and linking against an outside copy of all imported software that is externally developed and maintained |Full| ParMOO does not contain any other package's source code. Note that Python packages are imported using the conventional `sys.path` system. Alternative instances of a package can be used, for example, by including them through an appropriate definition of the `PYTHONPATH` environment variable.|
| 8. | Have no hardwired print or IO statements that cannot be turned off |Full| In ParMOO, errors are handled through Python's exception handling. The only hardwired print statements are in the examples and, optionally, through Python's info-level logging channel.|

### Recommended Policies

| # | Policy                 |Support| Notes                   |
|---|------------------------|-------|-------------------------|
|**R1.**| Have a public repository. |Full| https://github.com/parmoo/parmoo is publicly available. |
|**R2.**| Free all system resources acquired as soon as they are no longer needed. |Full| Python has built-in garbage collection that frees memory when it becomes unreferenced. |
|**R3.**| Provide a mechanism to export ordered list of library dependencies. |Full| The dependencies for ParMOO are given in `setup.py`. `pip install parmoo` installs required dependencies if they do not already exist. |
|**R4.**| Document versions of packages that it works with or depends upon, preferably in machine-readable form.  |Full| A full list of tested external dependencies is provided in [`REQUIREMENTS`](https://github.com/parmoo/parmoo/blob/main/REQUIREMENTS). |
|**R5.**| Have SUPPORT, LICENSE, and CHANGELOG files in top directory.  |Partial| All files are included in the main branch of the repository. |
|**R6.**| Have sufficient documentation to support use and further development. |Full| ParMOO provides documentation through a *Sphinx* framework. It is published on [readthedocs](https://parmoo.readthedocs.io), which includes a user guide covering quick-start, installation, and many usage details. There are tutorials and examples. The developer guide contains information on internal modules. |
|**R7.**| Be buildable using 64-bit pointers; 32-bit is optional |Full| There is no explicit use of pointers in ParMOO since Python handles pointers internally and depends on the install of Python, which will generally be 64-bit on supported systems.|
|**R8.**| Do not assume a full MPI communicator; allow for user-provided MPI communicator |Full| ParMOO can be run in parallel using [libEnsemble](https://github.com/Libensemble/libensemble). libEnsemble takes an MPI communicator as an option; see full details about the [MPIExecutor in libEnsemble](https://libensemble.readthedocs.io). |
|**R9.**| Use a limited and well-defined name space (e.g., symbol, macro, library, include) |Full| The `MOOP` class is the fundamental data structure in ParMOO; see the [diagram](https://parmoo.readthedocs.io/en/latest/_images/moop-uml.svg). ParMOO modules are located in the `parmoo` directory.|
|**R10.**| Give best effort at portability to key architectures |Full| ParMOO is being regularly tested on Unix/Linux and MacOS. The current set of automatically tested, common architectures is viewable [here](https://github.com/parmoo/parmoo/blob/main/.github/workflows/parmoo-ci.yml). |
|**R11.**| Install headers and libraries under `<prefix>/include` and `<prefix>/lib`, respectively |Full| The standard Python installation is used for Python dependencies. This installs external Python packages under `<install-prefix>/lib/python<X.Y>/site-packages/`.|
|**R12.**| All BAND compatibility changes should be sustainable |Full| The BAND-SDK-compatible package is in the standard release path. All the changes here should be sustainable.|
|**R13.**| Respect system resources and settings made by other previously called packages |Full| ParMOO does not modify system resources or settings.|
|**R14.**| Provide a comprehensive test suite for correctness of installation verification. |Full| In addition to the installation test noted above in M2, ParMOO contains a comprehensive set of unit tests that can be run, individually or all at once, via pytest with a high coverage. Running the provided ``.\run-tests.sh`` performs comprehensive testing.|
