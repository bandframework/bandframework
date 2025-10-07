'''
    Calculates the 12C(p,gamma) cross section and compares it to the Vogl data.
    "Free" parameters:
        * ANC (1/2-)
        * level energy (1/2+)
        * partial width (1/2+, elastic)
        * partial width (1/2+, capture)
'''

import os
import sys
from multiprocessing import Pool

import emcee
import numpy as np
from scipy import stats

import model

########################################
# We'll set up the sampler and get it started.

nd = model.azr.config.nd
nw = 2*nd # number of walkers = 2 * number of sampled parameters

# Pick a point (theta) in parameter space around which we'll start each walker.
theta0 = [2.1, 2.37, 33600, -0.6325, 1, 1]
# theta0 = np.array(model.azr.config.get_input_values())
# theta0 = np.array([[pi.rvs() for pi in model.priors] for _ in range(nw)])
# Each walkers needs its own starting position. We'll take normally distributed
# random values centered at theta0.
p0 = np.zeros((nw, nd))

mask = np.array([0.01, 0.0001, 0.01, 0.01, 0.01, 0.01])
for i in range(nw):
    mu = theta0
    sig = np.abs(theta0) * mask # 1% width
    p0[i, :] = stats.norm(mu, sig).rvs()

# We'll store the chain in test_mcmc.h5. (See emcee Backends documentation.)
backend = emcee.backends.HDFBackend('chain.h5')
backend.reset(nw, nd)

nsteps = 10 # How many steps should each walker take?
nthin = 1 # How often should the walker save a step?
nprocs = 4 # How many Python processes do you want to allocate?
# AZURE2 and emcee are both parallelized. We'll restrict AZURE2 to 1 thread to
# simplify things.
os.environ['OMP_NUM_THREADS'] = '1'

# emcee allows the user to specify the way the ensemble generates proposals.
moves = [(emcee.moves.DESnookerMove(), 0.8), (emcee.moves.DEMove(), 0.2)]

with Pool(processes=nprocs) as pool:
    sampler = emcee.EnsembleSampler(nw, nd, model.lnP, moves=moves, pool=pool,
            backend=backend)
    state = sampler.run_mcmc(p0, nsteps, thin_by=nthin, progress=True, tune=True)
