BAND Camp 2021: Introduction to Bfrescox
=========================================

:Authors: **Ozge Surer (Northwestern), Filomena Nunes (Michigan State), Matt Plumlee (Northwestern), Stefan Wild (Argonne/Northwestern)**
:Date: **12/13/2021**

This Bfrescox tutorial has four sections to illustrate how `frescox` interfaces with `surmise`.

In this tutorial, we use a Colab notebook, which allows you to run the code in an interactive and consistent environment. You can then modify this code and run it on your preferred system(s).


* `Tutorial I-Section I <BANDCamp_nbs/Bfrescox_intro.ipynb>`_ : To verify that the code works as expected.

* `Tutorial I-Section II <BANDCamp_nbs/Bfrescox_fit.ipynb>`_ : To illustrate the :math:`\\\chi^2` minimization wrapper for `frescox`.

* `Tutorial I-Section III <BANDCamp_nbs/Bfrescox_GPR.ipynb>`_ : A brief introduction to uncertainty quantification.

* `Tutorial I-Section IV <BANDCamp_nbs/Bfrescox_surmise.ipynb>`_ : To illustrate how `frescox` interfaces with `surmise`.

As of BAND Framework v0.5.0, the ptemcee sampler used in Section III is no longer compatible with the current numpy package.
We provide an updated notebook, `Tutorial I-Section III (updated) <BANDCamp_nbs/Bfrescox_GPR_v05.ipynb>`_ , which uses the emcee sampler instead.
