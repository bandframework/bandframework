# BAND SDK v0.2 Community Policy Compatibility for OpenBT

To suggest changes to these requirements or obtain more information, please contact [BAND](https://bandframework.github.io/team).

Details on citing the current version of the BAND Framework can be found in the [README](https://github.com/bandframework/bandframework).

**Website:** https://github.com/bandframework/OpenBT \
**Contact:** mpratola@iu.edu, jcyannotty@gmail.com \
**Icon:** https://github.com/bandframework/OpenBT/blob/main/docs/images/openbt_logo_rect.png \
**Description:**  OpenBT implements a variety of Bayesian tree models, including regression, model mixing, sensitivity analysis and multiobjective optimization.  It is built around a C++ code base along with Python and R interfaces.

### Mandatory Policies

**BAND SDK**
| # | Policy                 |Support| Notes                   |
|---|-----------------------|-------|-------------------------|
| 1. | Support BAND community GNU Autoconf, CMake, or other build options. |Full| The C++ package is built via Meson, while the Python package can be built and installed using pip.  The R package is installed from within R using the R CRAN "remote" package, which installs directly from the github repo. [M1 details](#m1-details) |
| 2. | Have a README file in the top directory that states a specific set of testing procedures for a user to verify the software was installed and run correctly. |Full| The README exists and clearly refers users to the User Guides for information on testing an installation. |
| 3. | Provide a documented, reliable way to contact the development team. |Full| The OpenBT team can be contacted via the listed author emails or via the public issues page on GitHub. |
| 4. | Come with an open-source license. |Full| Uses the MIT license. |
| 5. | Provide a runtime API to return the current version number of the software. |Partial| Users of the Python package can retrieve the current version using `openbt.__version__`. To be added in a future release of the C++ command line tools and R package. |
| 6. | Provide a BAND team-accessible repository. |Full| https://github.com/bandframework/OpenBT |
| 7. | Must allow installing, building, and linking against an outside copy of all imported software that is externally developed and maintained. |Full| The only external dependency for the C++ package is Eigen; the Meson build allows linking to a user-provided Eigen install or will install Eigen itself if no installation is detected during build. |
| 8. | Have no hardwired print or IO statements that cannot be turned off. |Partial| Building with Meson option `verbose=false` limits hardwired print statements. |

M1 details <a id="m1-details"></a>: The Meson build system allows one to build a stand-alone install of the base C++ package independent of Python or R.  In principle it can be used directly via the command line interface.  The Python package follows standard Python install procedures with pip pulling information from OpenBT's pypi.org entry.  By default the Python package also builds the base C++ package to facilitate simpler installation for Python users.  The R package requires independent installation of the C++ package first and then R users can directly install the Ropenbt interface from within R.


### Recommended Policies

| # | Policy                 |Support| Notes                   |
|---|------------------------|-------|-------------------------|
|**R1.**| Have a public repository. |Full| OpenBT is a public repository. |
|**R2.**| Free all system resources acquired as soon as they are no longer needed. |Full| None. |
|**R3.**| Provide a mechanism to export ordered list of library dependencies. |None| None. |
|**R4.**| Document versions of packages that it works with or depends upon, preferably in machine-readable form. |Partial| Python dependencies listed in `pyproject.toml`; R dependencies listed in `NAMESPACE`. |
|**R5.**| Have SUPPORT, LICENSE, and CHANGELOG files in top directory.  |Partial| LICENSE file provided. |
|**R6.**| Have sufficient documentation to support use and further development.  |Partial| A readthedocs.io site provides User and Developer Guides.  The Python and R User Guides provide examples. |
|**R7.**| Be buildable using 64-bit pointers; 32-bit is optional. |Partial| The C++ package is built for 64-bit. |
|**R8.**| Do not assume a full MPI communicator; allow for user-provided MPI communicator. |Full| While the C++ code can be built with MPI, it is compiled into standalone command line tools that each run MPI independently on each invocation.  The Python and R packages call these tools individually with the desired MPI parallelization setup.  The C++ code is not built into larger applications where fine-grained programmatic setting of MPI communicators is used. |
|**R9.**| Use a limited and well-defined name space (e.g., symbol, macro, library, include). |Full| None. |
|**R10.**| Give best effort at portability to key architectures. |Full| The Meson build allows installation on Linux and macOS and is tested regularly with Open MPI and MPICH.  Actions test the Python package similarly and across all supported versions of Python. |
|**R11.**| Install headers and libraries under `<prefix>/include` and `<prefix>/lib`, respectively. |Full| None. |
|**R12.**| All BAND compatibility changes should be sustainable. |Full| None. |
|**R13.**| Respect system resources and settings made by other previously called packages. |Full| None. |
|**R14.**| Provide a comprehensive test suite for correctness of installation verification. |Partial| Python package provides a test suite, which is used in test actions.  Some rudimentary testing of C++ package on build.  R package has no integrated testing due to R's package build/test runtime limitation.  To be improved in a future release. |
