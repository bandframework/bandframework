'''
Defines the Bayesian model we will use to analyze the Vogl and Rolfs&Azuma
12C(p,gamma) data.
'''

import numpy as np
from scipy import stats

from brick.azr import AZR

# With a preconfigured .azr file, we have all of the information we need to
# instantiate our AZR object.
azr = AZR('12C+p.azr')

# This tells the azr object where to "work". Instead of possibly cluttering up
# the local directory, azr will read/write to the /tmp director. (This also
# makes it easier to use with multi-node clusters.)
azr.root_directory = '/tmp/'

# We'll read the data from the output file since it's already in the
# center-of-mass frame.
data = np.loadtxt('output/AZUREOut_aa=1_R=2.out')
x = data[:, 0]  # energies
y_default_norm = data[:, 5]  # cross sections
dy_default_norm = data[:, 6]  # cross section uncertainties

########################################
# Next, let's set up the Bayesian calculation. Recall:
# * lnP \propto lnL + lnPi
# where
# * P = posterior
# * L = likelihood
# * Pi = prior

# We'll work from right to left.
# First, we need prior disributions for each sampled parameters.
priors = [
    stats.uniform(0, 5),
    stats.uniform(1, 5),
    stats.uniform(0, 50000),
    stats.uniform(-100, 200),
    stats.lognorm(0.1),
    stats.lognorm(0.1)
]

def lnPi(theta):
    '''
    Takes an input parameter array, theta, and returns the sum over the
    ln(probability density function) for each parameter based on the prior
    distributions defined in the "priors" list.
    '''
    return np.sum([pi.logpdf(t) for (pi, t) in zip(priors, theta)])


# To calculate the likelihood, we generate the prediction at theta and compare
# it to data. (Assumes data uncertainties are Gaussian and IID.)
def lnL(theta):
    '''
    Takes an input parameter array, theta, and returns the ln(Likelihood) based
    on the assumption that the data uncertainties are Gaussian and IID.
    '''
    output = azr.predict(theta)[0]
    mu = output.xs_com_fit
    y = output.xs_com_data
    dy = output.xs_err_com_data
    return np.sum(
        -np.log(np.sqrt(2*np.pi)*dy_default_norm) -
        0.5*((y - mu)/dy)**2
    )


def lnP(theta):
    '''
    The ln(Posterior) is simply the sum of the ln(Prior) and ln(Likelihood)
    functions.
    If any of the parameters fall outside of their prior distributions, go
    ahead and return lnPi = -infty. Don't bother running AZURE2 or risking
    calling it with a parameter value that will throw an error.
    '''
    lnpi = lnPi(theta)
    if lnpi == -np.inf:
        return lnpi
    return lnL(theta) + lnpi

