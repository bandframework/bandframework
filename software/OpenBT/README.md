# OpenBT in BAND 
OpenBT is an extensible software project managed in a dedicated [bandframework
repository](https://github.com/bandframework/OpenBT) that implements a variety
of Bayesian tree models for scientific and industry applications, including
regression, model mixing, sensitivity analysis and multiobjective optimization.
It also provides functionality to the Trees mixing interface of the Taweret
software package, which is also part of the BAND framework.

The heart of OpenBT is a set of C++ tools that can be used directly *via* the
command line or indirectly through the ``openbt`` Python package or ``RopenBT``
R package, which wrap the tools.  Typically these tools are built with an
implementation of the Message Passing Interface (MPI), such as Open MPI or
MPICH, to enable distributed parallelization of computations.

Please refer to our [User Guides](https://openbt.readthedocs.io) for more
information about installing, testing, and using the different OpenBT software
tools.

General details such as support, contributing, copyright, license, and citing OpenBT
are provided in the project's [landing
page](https://github.com/bandframework/OpenBT).

The [BAND SDK v0.2 Community Policy](/resources/sdkpolicies/bandsdk.md)
compatibility documentation for OpenBT is provided in
[OpenBTbandsdk.md](/software/OpenBT/OpenBTbandsdk.md).
