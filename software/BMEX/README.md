![](https://raw.githubusercontent.com/massexplorer/bmex-static/main/public/bmex-ico.png)

# The Bayesian Mass Explorer

The Bayesian Mass Explorer (BMEX) is a user-focused web application that provides a one-stop-shop for quantified theoretical model predictions of nuclear masses and related quantities.

A [BAND SDK v0.2 Community Policy](/resources/sdkpolicies/bandsdk.md) compatibility documentation for BMEX is contained in [BMEX-bandsdk.md](/software/BMEX/BMEX-bandsdk.imd).


## Accessing BMEX

While BMEX can be self-hosted and accessed locally, most users will want to access the tool through a modern web browser at [bmex.dev](https://bmex.dev).
From the landing page, links to the various sub-applications can be found, with the masses tool being the only public one at the moment.

One can cite the software with the following bibtex entry:
```
@software{bmex,
  author       = {Kyle Godbey, Landon Buskirk, and
                  Giuliani, Pablo},
  title        = {{BMEX} - {T}he {B}ayesian {M}ass {E}xplorer},
  month        = sep,
  year         = 2023,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.7111988},
  url          = {https://doi.org/10.5281/zenodo.7111988}
}
```

## Using BMEX

The intended usage of BMEX is to offer a user-friendly interface for the exploration of data and the extraction of information for nuclear theorists and experimentalists.

The most basic use case is to simply print the predictions of masses and derived from various theoretical models. This is done by selecting `Single Nucleus` in the `Dimension` drop-down on the left sidebar. You can then define your nucleus, quantity, and dataset to view what's in the BMEX database. If you select `All`, you will get a full dump of the mass related data available for that dataset.

While this is useful as a reference, there are also a suite of plotting tools for isotopic, isotonic, and isobaric chains (`1D Chains` in the drop-down) as well as entire nuclear landscapes.
The chain functionality is particularly useful for comparing how various quantities change as you change the chains of interest and for comparing multiple models and experimental data simultaneously.
The landscape feature is intended to provide a global view of these indicators for experimental and theoretical datasets to identify structure across the nuclear chart and, in the case of theory predictions beyond current experimental limits, highlight interesting regions that can be the subject of future experimental studies.

To manage the user's workspace, they can choose to either export the figures in the current view as vector PDFs or they can choose to save a link to the current view, allowing them to return at a future time or to share the workspace with others.

## Additional Support

For additional support, please use the issues tab of the [BMEX Github Repository](https://github.com/massexplorer/bmex-masses) to suggest features, report bugs, or get help with the tool. We also welcome all contributions via pull requests on the repository.
