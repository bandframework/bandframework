# frescox Installation

This document describes how to install `frescox` to run a scattering code for
coupled-channels calculations

`frescox` is publicly available at the repo https://github.com/LLNL/Frescox; more detailed documentation for earlier versions of `frescox` can be found at http://www.fresco.org.uk

In order to install `frescox`:

- Clone or download the repo  https://github.com/LLNL/Frescox to your machine

- Go to `Frescox/source` directory, in which you will find a `makefile` 

- Open `makefile` and remove one of comment prefix `#`s between lines 39-60 to select appropriate machine type

  - `MACH=i386` and `MACH=intel` have been BAND tested for macOS and intel machines, and accurately work for these architectures

- To compile `frescox`, go to `Frescox/source` and issue the following from your terminal:

```python
  make
  make install
  make clean
```
  - If you obtain an error running `make`:
  
    - In line 69, there is an `include fx$(MACH).def` to include a definition file in the source directory for the previously selected machine architecture; you may find that running the commands below produces an error such as 
`makefile:69: .def: No such file or directory`

      - If this `include` does not work for your machine, comment out `include` by adding a `#` prefix, and then copy all lines in `fx$(MACH).def` (for your specified `MACH` above) to the `makefile` starting from the line where `include` is commented out

      - As an example if you select `MACH=i386` on lines 39-60, then you should copy all lines in `fxi386.def` file to the `makefile`

    - In line 72 (of the original, unedited `makefile`), there is a line `LOCAL = f14$(MACH)`; you may find that you need to edit this line (e.g., because you receive a `fatal error: no input files; unwilling to write output files`)

      - If this line does not work with your architecture, replace `$(MACH)` in this line with your `MACH` defined above

      - As an example, if you have `MACH=i386`, then replace `LOCAL = f14$(MACH)` with `LOCAL = f14i386`

  - If you obtain an error `No such file or directory` running `make install`:
  
    - It assumed that the directory `$(HOME)/binw` exists; if it does not, you will need to create this directory (and probably the `MACH`-specific directory within it). Alternatively, you can replace references in `makefile` to `$(HOME)/binw` to point to whatever location in which you want the frescox executables to reside.

- To verify you have successfully compiled `frescox`, go to `Frescox/test` directory from your terminal and run one of the tests provided in that directory

  - As an example, type `frescox < lane20.nin > lane20.out` to run a test with input file `lane20.nin`, and the output file is written to `lane20.out`
  
  - If you obtain an error such as `frescox: command not found`, this means that the `frescox` executable is not in your path. You should add the location of this executable (see above reference to a subdirectory of `$(HOME)/binw`)

- In order to verify the code works as expected and perform Bayesian calibration, go to https://github.com/bandframework/privateband/tree/team/Software/Frescox/Tutorial_I
