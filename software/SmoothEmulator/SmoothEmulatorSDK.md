# BAND SDK v0.2 Community Policy Compatibility for SmoothEmulator


> This document summarizes the efforts of current and future BAND member packages to achieve compatibility with the BAND SDK community policies.  Additional details on the BAND SDK are available [here](https://github.com/bandframework/bandframework/blob/main/resources/sdkpolicies/bandsdk.md) and should be considered when filling out this form. The most recent copy of this template exists [here](https://github.com/bandframework/bandframework/blob/main/resources/sdkpolicies/template.md).
>
> This file should filled out and placed in the directory in the `bandframework` repository representing the software name appended by `bandsdk`.  For example, if you have a software `foo`, the compatibility file should be named `foobandsdk.md` and placed in the directory housing the software in the `bandframework` repository. No open source code can be included without this file.
>
> All code included in this repository will be open source.  If a piece of code does not contain a open-source LICENSE file as mentioned in the requirements below, then it will be automatically licensed as described in the LICENSE file in the root directory of the bandframework repository.
>
> Please provide information on your compatibility status for each mandatory policy and, if possible, also for recommended policies. If you are not compatible, state what is lacking and what are your plans on how to achieve compliance. For current BAND SDK packages: If you were not fully compatible at some point, please describe the steps you undertook to fulfill the policy. This information will be helpful for future BAND member packages.
>
> To suggest changes to these requirements or obtain more information, please contact [BAND](https://bandframework.github.io/team).


**Website:** https://github.com/bandframework/bandframework/blob/main/software/SmoothEmulator \
**Contact:** For assistance contact Scott Pratt (prattsc@msu.edu) \
**Icon:** N/a \
**Description:**  SmoothEmulator evaluates training points, fits a smooth emulator, and provides ways to explore a parameter space based on this emulator. 

### Mandatory Policies

**BAND SDK**
| # | Policy                 |Support| Notes                   |
|---|-----------------------|-------|-------------------------|
| 1. | Support BAND community GNU Autoconf, CMake, or other build options. |Full| All software is compiled through CMake. |
| 2. | Have a README file in the top directory that states a specific set of testing procedures for a user to verify the software was installed and run correctly. |Somewhat| User can run tutorial and check against expected output. |
| 3. | Provide a documented, reliable way to contact the development team. |Full| Yes, in README.|
| 4. | Come with an open-source license. |Full| Uses a [GPL v3 license](LICENSE.md).|
| 5. | Provide a runtime API to return the current version number of the software. |None| Unlikely to be implemented.|
| 6. | Provide a BAND team-accessible repository. |Full| Smooth Emulator is available via https://github.com/bandframework/bandframework. |
| 7. | Must allow installing, building, and linking against an outside copy of all imported software that is externally developed and maintained. |Full| No outside software.|
| 8. | Have no hardwired print or IO statements that cannot be turned off. |Full| Parameter files provide option to redirect output to file.|

### Recommended Policies

| # | Policy                 |Support| Notes                   |
|---|------------------------|-------|-------------------------|
|**R1.**| Have a public repository. |Partial| https://github.com/bandframework/bandframework is publicly available. |
|**R2.**| Free all system resources acquired as soon as they are no longer needed. |Not Checked| On todo list. |
|**R3.**| Provide a mechanism to export ordered list of library dependencies. |None| Not planned, beyond what CMake reports. |
|**R4.**| Document versions of packages that it works with or depends upon, preferably in machine-readable form.  |Full| Only package needed is Eigen, CMake sets version. |
|**R5.**| Have SUPPORT, LICENSE, and CHANGELOG files in top directory. |Partial| No CHANGELOG. |
|**R6.**| Have sufficient documentation to support use and further development. |Partial| Manual does not describe development. |
|**R7.**| Be buildable using 64-bit pointers; 32-bit is optional. |Full| Standard C++. |
|**R8.**| Do not assume a full MPI communicator; allow for user-provided MPI communicator. |N/a| Will consider parallelism for future. |
|**R9.**| Use a limited and well-defined name space (e.g., symbol, macro, library, include). |None| Not yet implemented.|
|**R10.**| Give best effort at portability to key architectures. |Full| Should work with any standard UNIX/Linux-based system (no Windows support provided).|
|**R11.**| Install headers and libraries under `<prefix>/include` and `<prefix>/lib`, respectively. |Full| All include paths have the form include/msu_*.|
|**R12.**| All BAND compatibility changes should be sustainable. |Full| The BAND-SDK-compatible package is in the standard release path. All the changes here should be sustainable.|
|**R13.**| Respect system resources and settings made by other previously called packages. |Not Implemented| May need to add some "delete"s.|
|**R14.**| Provide a comprehensive test suite for correctness of installation verification. |Partial| Only test is through running the provided tutorial.|
