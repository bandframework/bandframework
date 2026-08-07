Sunsetting support for elements of the BAND software framework
===============
The main purpose of this document is to provide guidance to BAND software authors who wish to cease maintaining their code. 

The sunsetting process has two main goals:
 1. To limit incorrect usage of the package by explicitly communicating that the package will no longer be maintained via all means that users might find the package
 2. To clearly communicate in official BAND framework documentation that the package to sunset was a part of the framework but is no longer maintained.

As an example of satisfying the first goal, consider a development team that no longer plans to maintain their Python package, which is publicly available on PyPI.  Some sensible steps to communicate responsibly and effectively with the user community are to
 1. Add text at the top of the root README in their package's git repository explaining that the software is no longer maintained and optionally providing suggestions/alternatives;
 2. Add similar remarks to the landing pages of public documentation sites such as Read the Docs, Jupyter books, and GitHub pages (or disable these);
 3. Update the README of their Python package to indicate the same;
 4. Update the internal short description in the package's definition to indicate that the package is deprecated, was moved, or was replaced;
 5. Change the status in their package's metadata (*e.g.,* ``Development Status :: 7 - Inactive``);
 6. Update the code so that it emits a deprecation warning when the package is imported;
 7. Document the last working versions of all external dependencies of the package (*e.g.,* use ``numpy>=1.0.1,<2.0.0`` instead of ``numpy>=1.0.1`` in their ``requirements.txt`` file); and
 8. Publish a final version of the package to PyPI so that the package's PyPI landing page is rendered with some of the above changes and any users that ``pip install`` the package without reading the documentation see the deprecation warning.

Hopefully, the steps needed to sunset non-Python projects can be derived from the previous example.  In such cases, please add final details here to help others with this process.

The second goal will be achieved as a natural part of a subsequent release of the BAND framework.  In particular,
as dictated by the `BAND Framework release process <https://github.com/bandframework/bandframework/blob/main/resources/dev_guide/release-proc.rst>`_, the team conducting the release will inventory all BAND components with the goal of
listing separately all packages that are no longer actively maintained.
