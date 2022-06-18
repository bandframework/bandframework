# BAND SDK v0.1 Community Policy Compatibility for Frescox


> This document summarizes the efforts of current and future BAND member packages to achieve compatibility with the BAND SDK community policies.  Additional details on the BAND SDK are available [here](/resources/sdkpolicies/bandsdk.md) and should be considered when filling out this form. The most recent copy of this template exists [here](/resources/sdkpolicies/template.md).
>
> This file should filled out and placed in the directory in the `bandframework` repository representing the software name appended by `bandsdk`.  For example, if you have a software `foo`, the compatibility file should be named `foobandsdk.md` and placed in the directory housing the software in the `bandframework` repository. No open source code can be included without this file.
>
> All code included in this repository will be open source.  If a piece of code does not contain a open-source LICENSE file as mentioned in the requirements below, then it will be automatically licensed as described in the LICENSE file in the root directory of the bandframework repository.
>
> Please provide information on your compatibility status for each mandatory policy and, if possible, also for recommended policies. If you are not compatible, state what is lacking and what are your plans on how to achieve compliance. For current BAND SDK packages: If you were not fully compatible at some point, please describe the steps you undertook to fulfill the policy. This information will be helpful for future BAND member packages.
>
> To suggest changes to these requirements or obtain more information, please contact [BAND](https://bandframework.github.io).


**Website:** https://github.com/LLNL/Frescox

**Contact:** nunes@frib.msu.edu

**Icon:** No icon info

**Description:** Scattering code Frescox for coupled-channels calculations 


### Mandatory Policies

**BAND SDK**

| # | Policy                 |Support| Notes                   |
|---|-----------------------|-------|-------------------------|
| 1. | Support BAND community GNU Autoconf, CMake, or other build options |Full|  |
| 2. | Provide a comprehensive test suite for correctness of installation verification |Full| The test/ directory contains : at least 6 test jobs xeta, lane20 & f19xfr, e80f49b, on2 & be11 their various outputs SUN/*.out |
| 3. | Provide a documented, reliable way to contact the development team |Full| |
| 4. | Come with an open-source license |Full| Frescox uses GNU general public license. |
| 5. | Provide a runtime API to return the current version number of the software |Full| |
| 6. | Provide a BAND team-accessible repository |Full| https://github.com/LLNL/Frescox |
| 7. | Must allow installing, building, and linking against an outside copy of all imported software that is externally developed and maintained |Full| |
| 8. |  Have no hardwired print or IO statements that cannot be turned off |None.| |



### Recommended Policies

| #  | Policy                 |Support| Notes                   |
|---|------------------------|-------|-------------------------|
|**R1.**| Have a public repository. |Full| https://github.com/LLNL/Frescox is publicly available. |
|**R2.**| Free all system resources acquired as soon as they are no longer needed. |Full|  |
|**R3.**| Provide a mechanism to export ordered list of library dependencies. |Full| |
|**R4.**| Document versions of packages that it works with or depends upon, preferably in machine-readable form.  |Full| |
|**R5.**| Have README, SUPPORT, LICENSE, and CHANGELOG files in top directory.  |Full| All files are included in the repository. |
|**R6.**| Have sufficient documentation to support use and further development. |Full| The man/ directory contains the instruction manual in latex: frescox-input-manual.tex: latex source frescox-input-manual.pdf: printable output More documentation is at http://www.fresco.org.uk/documentation.htm |
|**R7.**| Be buildable using 64-bit pointers; 32-bit is optional |Full| |
|**R8.**| Do not assume a full MPI communicator; allow for user-provided MPI communicator |N/a| None. |
|**R9.**| Use a limited and well-defined name space (e.g., symbol, macro, library, include) |Full| |
|**R10.**| Give best effort at portability to key architectures |Full|  |
|**R11.**| Install headers and libraries under `<prefix>/include` and `<prefix>/lib`, respectively |Full| |
|**R12.**| All BAND compatibility changes should be sustainable |Full| |
|**R13.**| Respect system resources and settings made by other previously called packages |Full| |
