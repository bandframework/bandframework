Sunsetting elements of the BAND software Framework
===============
The main purpose of this document is to provide guidance to BAND software authors who wish to cease maintaining their code. 

If the development team no longer plans to maintain the package: 
 1. Text must be added to top-level repository README.md explaining this;
 2. The last working versions of the dependencies of package should be documented (e.g. by using `numpy>=1.0.1,<2.0.0` instead of `numpy>=1.0.1` in the requirements.txt). 

As part of the BAND Framework release process, the team conducting the release will then inventory all BAND components, with the goal of
listing the ones that are not being actively maintained, see bandframework/resources/release-proc.rst.
