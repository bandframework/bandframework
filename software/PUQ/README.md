# PUQ in BAND 

PUQ is a Python package for generating experimental designs tailored for uncertainty quantification and featuring parallel implementations.

A [BAND SDK v0.2 Community Policy](/resources/sdkpolicies/bandsdk.md) compatibility documentation for PUQ is contained in [PUQ-bandsdk.md](/software/PUQ/PUQ-bandsdk.md).

Complete installation and testing details for PUQ are available at the [PUQ repo](https://github.com/parallelUQ/PUQ). PUQ documentation is available on [readthedocs](https://puq.readthedocs.io/).

The easiest way to get an installation of PUQ is to follow the guidance at the [PUQ README](https://github.com/parallelUQ/PUQ/blob/main/README.rst).

## PUQ Tutorial: Sequential Design for the Calibration of 48Ca(n,n)48Ca Reaction

We demonstrate our sequential strategy using a nuclear physics model to predict differential cross sections as a function of angle. The angular cross sections vary with different parametrizations of the optical potential, which are inputs to the reaction code ``frescox``. ``frescox`` generates cross sections for angles ranging from 0° to 180°. This case study focuses on elastic scattering data from the 48Ca(n,n)48Ca  reaction to find the optimal parametrization of the optical potential. We illustrate how ``PUQ`` and ``frescox`` work together through an example.

Further details on this problem, including a benchmark comparing various acquisition functions and parallel implementations, can be found in Section 8.3 of the paper by [Sürer, Plumlee, and Wild, 2024](https://www.tandfonline.com/doi/abs/10.1080/00401706.2023.2246157?src=&journalCode=utch20).

After installing ``PUQ``, ``frescox`` must also be installed to collect data using the sequential procedure. Additional notes on obtaining and building ``frescox`` can be found in the [Bfrescox README](/software/Bfrescox/README.md).

From the root directory of ``PUQ``, navigate to the ``examples/fresco_example`` directory:
```
cd examples/fresco_example
```

The template input file ``48Ca_template.in`` is provided to run the example.

In this illustration, we acquire 64 data points (``-max_eval 64``), including an initial design of size 32 (``-n_init_thetas 32``) from a uniform prior. To acquire data points for calibrating the model via the expected integrated variance criterion (``al_func "eivar"``), run the following command:
```
python3 puq_bfresco_test.py -max_eval 64 -al_func "eivar" -n_init_thetas 32
```

Running this script takes approximately 150 seconds on a personal Mac laptop. Upon completion, the file ``Figure_fresco.jpg`` will be saved in the ``examples/fresco_example`` directory.

After collecting data, we observe the cross sections evaluated with the acquired parameters.

![Illustration of PUQ with Bfrescox](https://github.com/parallelUQ/PUQ/blob/main/examples/fresco_example/Figure_fresco.png)


