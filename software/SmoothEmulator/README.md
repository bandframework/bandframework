# Smooth Emulator

Smooth Emulator is a software project managed in a dedicated [bandframework
repository](https://github.com/bandframework/SmoothEmulator) that can be used to
emulate full models that are *smooth*.  Given a
* model to emulate,
* its list of model parameters,
* the relative importance of the parameters, and
* a description of their priors,

Smooth Emulator provides an optimal set of training points in the model's space.
The software constructs and tunes emulators for each observable of interest
using the model values and uncertainties derived from the training set.  The
project provides a Markov chain Monte Carlo sampler that can be used with the
emulators and combined with experimental results and uncertainties to acquire a
set of samples in the model's parameter space that are consistent with the
posterior distribution.

While Smooth Emulator is implemented in C++, a Python binding is also provided.

Please refer to the [User
Manual](https://github.com/bandframework/SmoothEmulator/blob/main/doc/UserManual.pdf)
for detailed descriptions and installation directions.  The User Manual also
provides a tutorial.

The [BAND SDK v0.2 Community Policy](/resources/sdkpolicies/bandsdk.md)
compatibility documentation for Smooth Emulator is provided in
[SmoothEmulatorSDK.md](/software/SmoothEmulator/SmoothEmulatorSDK.md).
