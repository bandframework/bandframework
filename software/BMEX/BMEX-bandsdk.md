# BAND SDK v0.2 Community Policy Compatibility for BMEX


> This document summarizes the efforts of current and future BAND member packages to achieve compatibility with the BAND SDK community policies.  Additional details on the BAND SDK are available [here](/resources/sdkpolicies/bandsdk.md) and should be considered when filling out this form. The most recent copy of this template exists [here](/resources/sdkpolicies/template.md).
>
> This file should filled out and placed in the directory in the `bandframework` repository representing the software name appended by `bandsdk`.  For example, if you have a software `foo`, the compatibility file should be named `foobandsdk.md` and placed in the directory housing the software in the `bandframework` repository. No open source code can be included without this file.
>
> All code included in this repository will be open source.  If a piece of code does not contain a open-source LICENSE file as mentioned in the requirements below, then it will be automatically licensed as described in the LICENSE file in the root directory of the bandframework repository.
>
> Please provide information on your compatibility status for each mandatory policy and, if possible, also for recommended policies. If you are not compatible, state what is lacking and what are your plans on how to achieve compliance. For current BAND SDK packages: If you were not fully compatible at some point, please describe the steps you undertook to fulfill the policy. This information will be helpful for future BAND member packages.
>
> To suggest changes to these requirements or obtain more information, please contact [BAND](https://bandframework.github.io/team).
>
> Details on citing the current version of the BAND Framework can be found in the [README](https://github.com/bandframework/bandframework).

**Website:** https://bmex.dev \
**Contact:** kyle@bmex.dev \
**Icon:** https://raw.githubusercontent.com/massexplorer/bmex-static/main/public/bmex-ico.png \
**Description:** The Bayesian Mass Explorer is a cloud-hosted web application that aims to provide a user-friendly interface to experimental data and theoretical model predictions of nuclear masses and related quantities.


### Mandatory Policies

**BAND SDK**
| # | Policy                 |Support| Notes                   |
|---|-----------------------|-------|-------------------------|
| 1. | Support BAND community GNU Autoconf, CMake, or other build options. |Full| BMEX is automatically built and deployed as a Docker container. When ran locally, BMEX is not required to be built and is entirely interpreted. |
| 2. | Have a README file in the top directory that states a specific set of testing procedures for a user to verify the software was installed and run correctly. | Full| The frontend will appear if the files and database have been downloaded appropriately. |
| 3. | Provide a documented, reliable way to contact the development team. |Full| The BMEX team can be contacted through the public [issues page on GitHub](https://github.com/massexplorer/bmex-masses/issues) or via an e-mail to [kyle@bmex.dev](kyle@bmex.dev).|
| 4. | Come with an open-source license |Full| The BMEX web application is released under the MIT license. |
| 5. | Provide a runtime API to return the current version number of the software. |Full| The BMEX web application is not interfaced in this way; there is no API endpoint at this time.|
| 6. | Provide a BAND team-accessible repository. |Full| https://github.com/massexplorer/bmex-masses |
| 7. | Must allow installing, building, and linking against an outside copy of all imported software that is externally developed and maintained |Full| BMEX does not contain other software's source code; any external dependencies are installed in the automated build process.|
| 8. | Have no hardwired print or IO statements that cannot be turned off |Full| For production deployments BMEX will only report fatal errors in the container log. For development, the web framework will report errors in the debug tools. |

### Recommended Policies

| # | Policy                 |Support| Notes                   |
|---|------------------------|-------|-------------------------|
|**R1.**| Have a public repository. |Full| https://github.com/massexplorer/bmex-masses is publicly available. |
|**R2.**| Free all system resources acquired as soon as they are no longer needed. |Full| Python has built-in garbage collection that frees memory when it becomes unreferenced. The user can choose to clear their browser cache after a BMEX session. |
|**R3.**| Provide a mechanism to export ordered list of library dependencies. |Full| The dependencies for BMEX are listed in requirements.txt and environment.yml, and are automatically compiled when the container is built.|
|**R4.**| Document versions of packages that it works with or depends upon, preferably in machine-readable form.  |Full| The dependencies for BMEX are listed in requirements.txt and environment.yml, and are automatically compiled when the container is built. |
|**R5.**| Have SUPPORT, LICENSE, and CHANGELOG files in top directory.  |Partial| LICENSE is in the root of the repository, a changelog is provided via git history, and support is encouraged through github issues. |
|**R6.**| Have sufficient documentation to support use and further development. |Partial| Code is documented, though no docs website currently exists. This is planned for future releases. |
|**R7.**| Be buildable using 64-bit pointers; 32-bit is optional. |Full| There is no explicit use of pointers in BMEX since Python handles pointers internally and depends on the install of Python, which will generally be 64-bit on supported systems.|
|**R8.**| Do not assume a full MPI communicator; allow for user-provided MPI communicator. |N/A| BMEX is a no MPI zone. |
|**R9.**| Use a limited and well-defined name space (e.g., symbol, macro, library, include) |Full| The BMEX web application is a self-contained application and, in production deployments, entirely separated from system installations of existing software.|
|**R10.**| Give best effort at portability to key architectures |Full| The flexible Plotly Dash web framework is portable thanks to its python backend. |
|**R11.**| Install headers and libraries under `<prefix>/include` and `<prefix>/lib`, respectively. |N/A| N/A |
|**R12.**| All BAND compatibility changes should be sustainable. |N/A| N/A|
|**R13.**| Respect system resources and settings made by other previously called packages. |Full| BMEX does not modify system resources or settings.|
|**R14.**| Provide a comprehensive test suite for correctness of installation verification. |N/A| N/A|
