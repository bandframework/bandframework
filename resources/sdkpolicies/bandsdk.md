## BAND Software Development Kit (SDK)

The BAND collaboration will carry out the software development in several stages that will include parallel lines of development and testing. External code delivery will happen via a [public Github repository](https://github.com/bandframework/bandframework). To suggest changes to these requirements or obtain more information, please contact members of [BAND](https://bandframework.github.io/team).

## Including code in BAND

All code included in the public Github `bandframework` repository will be open source.  If a piece of code does not contain a open-source LICENSE file as mentioned in the requirements below, then it will be automatically licensed as described in the LICENSE file in the root directory of the bandframework repository.  

New code can be included in the `bandframework` repository via a pull request to the `development` branch.  The name of the software should be the subdirectory name in the `/software` directory.  

BAND packages should include a compatibility document and a [template](https://raw.githubusercontent.com/bandframework/bandframework/main/resources/sdkpolicies/template.md) is provided. The compatibility  file should be placed in the root directory of the new subdirectory and labeled by the software name appended by SDKcompatibility.  For example, if you have a software `foo`, you should create a directory `/software/foo` and place the compatibility file named `fooSDKcompatibility.md` in the directory `/software/foo` alongside contributed code.  All pull requests for inclusion of a new `/software` directory  that do not include this document will be rejected.

## BAND SDK v0.1 Community Policies

The BAND Framework will conform with BAND Software Development Kit (SDK) requirements, summarized below. By using such requirements we envision ready interoperability across the BANDIT software ecosystem, large-scale scientific simulation codes, and other numerical libraries.

### Mandatory Policies

| # | Policy                 |
|---|-----------------------|
| 1. | Support BAND community GNU Autoconf, CMake, or other build options.
| 2. | Provide a comprehensive test suite for correctness of installation verification.
| 3. | Provide a documented, reliable way to contact the development team.
| 4. | Come with an open-source license.
| 5. | Provide a runtime API to return the current version number of the software.
| 6. | Provide a BAND team-accessible repository.
| 7. | Must allow installing, building, and linking against an outside copy of all imported software that is externally developed and maintained.
| 8. |  Have no hardwired print or IO statements that cannot be turned off.

### Recommended Policies

| # | Policy                 |
|---|------------------------|
|**R1.**| Have a public repository.
|**R2.**| Free all system resources acquired as soon as they are no longer needed.
|**R3.**| Provide a mechanism to export ordered list of library dependencies.
|**R4.**| Document versions of packages that it works with or depends upon, preferably in machine-readable form.
|**R5.**| Have README, SUPPORT, LICENSE, and CHANGELOG files in top directory.
|**R6.**| Have sufficient documentation to support use and further development.
|**R7.**| Be buildable using 64-bit pointers; 32-bit is optional.
|**R8.**| Do not assume a full MPI communicator; allow for user-provided MPI communicator.
|**R9.**| Use a limited and well-defined name space (e.g., symbol, macro, library, include).
|**R10.**| Give best effort at portability to key architectures.
|**R11.**| Install headers and libraries under `<prefix>/include` and `<prefix>/lib`, respectively.
|**R12.**| All BAND compatibility changes should be sustainable.
|**R13.**| Respect system resources and settings made by other previously called packages.
