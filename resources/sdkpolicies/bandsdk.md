## BAND Software Development Kit (SDK)

The BAND collaboration will carry out the software development in several stages that will include parallel lines of development and testing. External code delivery will happen via a [public Github repository](https://github.com/bandframework/bandframework). To suggest changes to these requirements or obtain more information, please contact members of [BAND](https://bandframework.github.io/team).

## Including code in BAND

All code included in the public Github `bandframework` repository will be open source.  If a piece of code does not contain a open-source LICENSE file as mentioned in the requirements below, then it will be automatically licensed as described in the [LICENSE](/LICENSE) file in the root directory of the `bandframework` repository.  

New code can be included in the `bandframework` repository via a pull request to the `develop` branch.  The name of the software should be the subdirectory name in the `/software` directory.  

BAND packages should include a compatibility document and a [template](/resources/sdkpolicies/template.md) is provided. The compatibility  file should be placed in the root directory of the new subdirectory and labeled by the software name appended by `bandsdk`.  For example, if you have a software `foo`, you should create a directory `/software/foo` and place the compatibility file named `foobandsdk.md` in the directory `/software/foo` alongside contributed code.  All pull requests for inclusion of a new `/software` directory  that do not include this document will be rejected.  You can also include a repository as a submodule in the `bandframework` repository by placing a pull request.  Please see the [additional instructions for including a submodule](/resources/dev_guide/git_instructions_for_submodules.md); all submodules must meet the same requirements as outlined here.

## BAND SDK v0.2 Community Policies

The BAND Framework will conform with BAND Software Development Kit (SDK) requirements, summarized below. By using such requirements we envision ready interoperability across the BAND software ecosystem, large-scale scientific simulation codes, and other numerical libraries.

As shown in the [template](/resources/sdkpolicies/template.md) and [linked SDK compatibility examples](/resources/sdkpolicies/README.md), we note that there are different levels to which packages can satisfy BAND's mandatory and recommended policies in order to facilitate interoperability with other BAND compatible tools. We emphasize that the most important thing is for packages to state their compliance with each policy.

### Mandatory Policies

| # | Policy                 |
|---|-----------------------|
| 1. | Support BAND community GNU Autoconf, CMake, or other build options.
| 2. | Have a README file in the top directory that states a specific set of testing procedures for a user to verify the software was installed and run correctly.
| 3. | Provide a documented, reliable way to contact the development team.
| 4. | Come with an open-source license.
| 5. | Provide a runtime API to return the current version number of the software.
| 6. | Provide a BAND team-accessible repository.
| 7. | Must allow installing, building, and linking against an outside copy of all imported software that is externally developed and maintained.
| 8. | Have no hardwired print or IO statements that cannot be turned off.

### Recommended Policies

| # | Policy                 |
|---|------------------------|
|**R1.**| Have a public repository.
|**R2.**| Free all system resources acquired as soon as they are no longer needed.
|**R3.**| Provide a mechanism to export ordered list of library dependencies.
|**R4.**| Document versions of packages that it works with or depends upon, preferably in machine-readable form.
|**R5.**| Have SUPPORT, LICENSE, and CHANGELOG files in top directory.
|**R6.**| Have sufficient documentation to support use and further development.
|**R7.**| Be buildable using 64-bit pointers; 32-bit is optional.
|**R8.**| Do not assume a full MPI communicator; allow for user-provided MPI communicator.
|**R9.**| Use a limited and well-defined name space (e.g., symbol, macro, library, include).
|**R10.**| Give best effort at portability to key architectures.
|**R11.**| Install headers and libraries under `<prefix>/include` and `<prefix>/lib`, respectively.
|**R12.**| All BAND compatibility changes should be sustainable.
|**R13.**| Respect system resources and settings made by other previously called packages.
|**R14.**| Provide a comprehensive, automated test suite for correctness of installation verification.

## Citing the BAND Framework

Please use the following to cite the BAND Framework:

    @techreport{bandframework,
        title       = {{BANDFramework: An} Open-Source Framework for {B}ayesian Analysis of Nuclear Dynamics},
        author      = {Kyle Beyer and Landon Buskirk and Moses Y-H. Chan and Tyler H. Chang and Richard James DeBoer and 
        Richard J. Furnstahl and Pablo Giuliani and Kyle Godbey and Kevin Ingles and Dananjaya Liyanage and Filomena M. Nunes and 
        Daniel Odell and Daniel R. Phillips and Matthew Plumlee and Matthew T. Pratola and Alexandra C. Semposki and \"Ozge S\"urer and 
        Stefan M. Wild and John C. Yannotty},
        institution = {},
        number      = {Version 0.3.0},
        year        = {2023},
        url         = {https://github.com/bandframework/bandframework}
    }
