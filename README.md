![BAND](resources/BAND_logo_v2.png)

# BAND Framework
This contains the primary public repository for the [BAND framework project](https://bandframework.github.io/). 

## Goals

The Bayesian Analysis for Nuclear Dynamics (BAND) software Framework provides tools and examples that 
facilitate principled Uncertainty Quantification in Nuclear Physics. 

This framework is funded by the NSF Office of Advanced Cyberinfrastructure, Cyberinfrastructure for Sustained Scientific Innovation program, under grant OAC-2004601.

We provide tools and examples that demonstrate how
- emulation of computationally expensive models,
- model calibration, and
- Bayesian model mixing

can be combined in order to provide a full accounting of the uncertainties in Nuclear Physics models–-including model
uncertainty.

These tools are designed to:
- enable predictions with full UQ for experimentally inaccessible environments such as neutron stars.
- when combined with code that implements Bayesian experimental design methodology, permit a quantitative assessment of the impact of proposed experiments.

More about the goals and structure of the BAND project can be found in D.R. Phillips et al., "[Get on the BAND Wagon: a Bayesian framework for quantifying model uncertainties in nuclear dynamics"](https://doi.org/10.1088/1361-6471/abf1df), J. Phys. G **48** (2021) 7, 072001.

A full list of BAND members together with a current list of  publications produced by the members using the tools and ideas of the project is available at https://bandframework.github.io


## BAND Framework Elements

BAND Framework elements are of two main types:
- BAND tools: these are pieces of python code that can be invoked to perform specific emulation, calibration, model-mixing, experimental-design, or linkage functions.
- BAND examples: these are typically notebooks, that show how the BAND tools can be used, singly or in combination, to quantify uncertainty in a specific nuclear-physics problem. 

BAND Framework tools and examples are found in [software/](/software/).

As of version 0.2.0+dev, the following tools are included:

- surmise: A surrogate model interface for calibration, uncertainty quantification, and sensitivity analysis.
- SaMBA: The Sandbox for Mixing via Bayesian Analysis.
- parMOO: A Python library for parallel multiobjective simulation optimization.
- Taweret: A Python package containing multiple Bayesian Model Mixing methods.

The following examples of the use of surmise are part of version 0.2.0+dev

- Bfrescox: A BAND extension of the frescox scattering code for coupled-channels calculations.

Version 0.2.0+dev also includes two examples of the use of Bayesian parameter estimation in nuclear-physics contexts:

- BRICK: the Bayesian R-matrix Inference Code Kit, facilitates extraction of R-matrix parameters from experimental data.
- QGP_Bayes: provides a tutorial on the use of JETSCAPE_SIMS tools to infer parameters of the QGP. 

## Downloading and using the BAND Framework

You are free to use any pieces of the BAND Framework that will advance your own research. Please cite the framework, and the original BAND paper, as detailed below under "Citing the BAND Framework".

## Submodules

BAND Framework currently includes some dependencies via git submodules. Currently, the following submodules are employed:

* [software/](software/)QGP_Bayes
* [software/](software/)SAMBA 
* [software/](software/)surmise
* [software/](software/)Taweret

As a consequence, when cloning the BAND Framework repository, the submodules can be retrieved automatically via
- `git clone --recursive` (in place of the usual `git clone`)

If you have already cloned the repository, the following modified `git` commands can be used:
- `git submodule update --init --recursive` (to obtain all submodules and any submodules those submodules have)
- `git pull --recurse-submodules=yes` (in place of the usual `git pull`)
- `git submodule update --init` (additionally, after `git pull`)
- `git submodule update --init software/surmise`
  (variant of the previous item, in case you want to only get the surmise submodule)

Note that submodules work modern git, version >= 2.38.0.

## Contributing to the BAND Framework

BAND welcomes contributions to the BAND Framework in a variety of forms; please see [CONTRIBUTING](CONTRIBUTING.rst).

The BAND Framework maintains a [BAND Software Development Kit (SDK)](/resources/sdkpolicies/bandsdk.md) that includes requirements and recommendations for contributing a package to the BAND Framework. 

Detailed instructions for contributing software, including by submodules, are found in BAND's [developer's guide](/resources/dev_guide).

## License 

All code included in the BAND Framework is open source, with the particular form of license contained in the top -level subdirectories of [software/](/software/).  If such a subdirectory does not contain a LICENSE file, then it is automatically licensed as described in the otherwise encompassing BAND Framework [LICENSE](/LICENSE).  

## Citing the BAND Framework

Please use the following to cite the BAND Framework:

    @techreport{bandframework,
        title       = {{BANDFramework: An} Open-Source Framework for {B}ayesian Analysis of Nuclear Dynamics},
        author      = {Moses Y-H. Chan and Richard James DeBoer and Richard J. Furnstahl and Dananjaya Liyanage and Filomena M. Nunes and 
        Daniel Odell and Daniel R. Phillips and Matthew Plumlee and Alexandra C. Semposki and \"Ozge S\"urer and Stefan M. Wild},
        institution = {},
        number      = {Version 0.2.0},
        year        = {2022},
        url         = {https://github.com/bandframework/bandframework}
    }
    
If possible, please also cite the original BAND Framework paper:

    @article{Phillips:2020dmw,
        author = "Phillips, D. R. and others",
        title = "{Get on the BAND Wagon: A Bayesian Framework for Quantifying Model Uncertainties in Nuclear Dynamics}",
        eprint = "2012.07704",
        archivePrefix = "arXiv",
        primaryClass = "nucl-th",
        doi = "10.1088/1361-6471/abf1df",
        journal = "J. Phys. G",
        volume = "48",
        number = "7",
        pages = "072001",
        year = "2021"
    }


## Resources
For more information, please see the [BAND framework website](https://bandframework.github.io/). 

Our [release process](resources/dev_guide/release-proc.rst) is also provided.
