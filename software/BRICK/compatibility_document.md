# BAND SDK v0.1 Community Policy Compatibility for BRICK


> This document summarizes the efforts of current and future BAND member
> packages to achieve compatibility with the BAND SDK community policies.
> Additional details on the BAND SDK are available
> [here](https://raw.githubusercontent.com/bandframework/bandframework/main/resources/sdkpolicies/bandsdk.md)
> and should be considered when filling out this form. The most recent copy of
> this template exists
> [here](https://raw.githubusercontent.com/bandframework/bandframework/main/resources/sdkpolicies/template.md).
>
> This file should filled out and placed in the directory in the `bandframework`
> repository representing the software name appended by `bandsdk`.  For example,
> if you have a software `foo`, the compatibility file should be named
> `foobandsdk.md` and placed in the directory housing the software in the
> `bandframework` repository. No open source code can be included without this
> file.
>
> All code included in this repository will be open source.  If a piece of code
> does not contain a open-source LICENSE file as mentioned in the requirements
> below, then it will be automatically licensed as described in the LICENSE file
> in the root directory of the bandframework repository.
>
> Please provide information on your compatibility status for each mandatory
> policy and, if possible, also for recommended policies. If you are not
> compatible, state what is lacking and what are your plans on how to achieve
> compliance. For current BAND SDK packages: If you were not fully compatible at
> some point, please describe the steps you undertook to fulfill the policy.
> This information will be helpful for future BAND member packages.
>
> To suggest changes to these requirements or obtain more information, please
> contact [BAND](https://bandframework.github.io/team).


**Website:** https://github.com/odell/brick
**Contact:** dodell@ohio.edu
**Icon:** 🧱
**Description:**  BRICK is a Python layer to the $R$-matrix code AZURE2. It
allows the user to easily sample what it typically optimized in an $R$-matrix
calculation. What is included here is a basic tutorial to familiarize BRICK
users with its most important features.

### Mandatory Policies

**BAND SDK**
| # | Policy                 |Support| Notes                   |
|---|-----------------------|-------|-------------------------|
| 1. | Support BAND community GNU Autoconf, CMake, or other build options. | Full | BRICK relies on the PyPA `build` package. |
| 2. | Provide a comprehensive test suite for correctness of installation verification. | Full | Basic functionality tests can be found in `tests/12Cpg` |
| 3. | Provide a documented, reliable way to contact the development team. |Full| In addition to the contact information above, BRICK is hosted on GitHub where the development team welcomes feedback and contributions from the community through issues and pull requests, respectively. |
| 4. | Come with an open-source license |Full| Uses MIT License. |
| 5. | Provide a runtime API to return the current version number of the software. |Full| The `brick` module has a `__version__` attribute that returns the version number as a string. |
| 6. | Provide a BAND team-accessible repository. |Full| [https://github.com/odell/brick](https://github.com/odell/brick) |
| 7. | Must allow installing, building, and linking against an outside copy of all imported software that is externally developed and maintained. | Full | |
| 8. | Have no hardwired print or IO statements that cannot be turned off. |Full| BRICK's main class, `AZR`, has a `verbose` attribute that can be set to `False` to silence all print statements. |

### Recommended Policies

| # | Policy                 |Support| Notes                   |
|---|------------------------|-------|-------------------------|
|**R1.**| Have a public repository. | Full | [https://github.com/odell/brick](https://github.com/odell/brick) |
|**R2.**| Free all system resources acquired as soon as they are no longer needed. | Full | |
|**R3.**| Provide a mechanism to export ordered list of library dependencies. | None |  |
|**R4.**| Document versions of packages that it works with or depends upon, preferably in machine-readable form.  | None |  |
|**R5.**| Have README, SUPPORT, LICENSE, and CHANGELOG files in top directory.  | Partial | Missing SUPPORT and CHANGELOG. |
|**R6.**| Have sufficient documentation to support use and further development.  | Partial |  |
|**R7.**| Be buildable using 64-bit pointers; 32-bit is optional. | Unsure | The PyPA `build` package may take care of this. |
|**R8.**| Do not assume a full MPI communicator; allow for user-provided MPI communicator. | N/A |  |
|**R9.**| Use a limited and well-defined name space (e.g., symbol, macro, library, include). | Partial | |
|**R10.**| Give best effort at portability to key architectures. | Full | |
|**R11.**| Install headers and libraries under `<prefix>/include` and `<prefix>/lib`, respectively. | N/A | |
|**R12.**| All BAND compatibility changes should be sustainable. | Full | |
|**R13.**| Respect system resources and settings made by other previously called packages. | Full | |
