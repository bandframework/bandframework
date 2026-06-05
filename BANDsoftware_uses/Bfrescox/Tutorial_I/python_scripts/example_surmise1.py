from frescox_wrapper import generate_input_file
from frescox_wrapper import frescox_output
import scipy.stats as sps
import numpy as np
import matplotlib.pyplot as plt
from surmise.emulation import emulator
from surmise.calibration import calibrator


def visualize(cal_object, emu_object, x_tr, theta_tr):
    import pylab

    theta_cal = cal_object.theta.rnd(1000)
    pred_cal = emu_object.predict(x=x_tr, theta=theta_cal)
    pred_cal_mean = pred_cal.mean()

    upper = np.percentile(pred_cal_mean, 97.5, axis=1)
    lower = np.percentile(pred_cal_mean, 2.5, axis=1)
    median = np.percentile(pred_cal_mean, 50, axis=1)

    titles = ['V', 'r', 'a', 'Ws', 'rs', 'as']

    fig, axes = plt.subplots(2, 3, figsize=(10, 6))
    k = 0
    for i in range(2):
        for j in range(3):
            axes[i, j].hist(theta_tr[:, k], alpha=0.5, color='r', label='prior')
            axes[i, j].hist(theta_cal[:, k], alpha=0.5, label='posterior')
            axes[i, j].set_title(titles[k])
            k += 1

    plt.legend(bbox_to_anchor=(-0.5, -0.25), loc='center', borderaxespad=0.1)
    plt.show()

    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(median, color='red', label='Prediction Mean')
    ax.plot(data[:, 0], data[:, 1], 'ko', label='Elastic Scattering Data')
    plt.xlabel(r'$\theta$ (deg)')
    plt.ylabel(r'd$\sigma/d\Omega$')
    ax.fill_between(range(0, 181), lower, upper, color='grey', alpha=0.5)
    plt.legend()
    ax.set_yscale('log')
    pylab.show()

# Define a class for prior of 6 parameters
class prior_fresco_Ca:
    """ This defines the class instance of priors provided to the method. """

    # 1.1, 0.05
    # Parameter names:  V, r, a, Ws, rs, as
    def lpdf(theta):
        return (sps.norm.logpdf(theta[:, 0], 49, 4.9) +
                sps.norm.logpdf(theta[:, 1], 0.90, 0.09) +
                sps.norm.logpdf(theta[:, 2], 0.67, 0.06) +
                sps.norm.logpdf(theta[:, 3], 3.39, 0.33) +
                sps.norm.logpdf(theta[:, 4], 1.08, 0.1) +
                sps.norm.logpdf(theta[:, 5], 0.27, 0.03)).reshape((len(theta), 1))

    def rnd(n):
        return np.vstack((sps.norm.rvs(49, 4.9, size=n),
                          sps.norm.rvs(0.90, 0.09, size=n),
                          sps.norm.rvs(0.67, 0.06, size=n),
                          sps.norm.rvs(3.39, 0.33, size=n),
                          sps.norm.rvs(1.08, 0.1, size=n),
                          sps.norm.rvs(0.27, 0.03, size=n))).T


def frescox_run(theta):
    theta_list = theta.tolist()
    f_l = []
    for para_obs in theta_list:
        generate_input_file(para_obs)
        exp_output = frescox_output()
        f_l.append(exp_output)

    return np.asarray(f_l).T


# Define x (input) and f (function evaluation)
theta_tr = prior_fresco_Ca.rnd(500)
f_tr = frescox_run(theta_tr)
x_tr = np.array(range(181)).reshape((181, 1))

# (Test) draw 50 random parameters from the prior
theta_test = prior_fresco_Ca.rnd(50)
f_test = frescox_run(theta_test)
x_test = np.array(range(181)).reshape((181, 1))

# Build an emulator for the frescox simulation
emu_fresco = emulator(x=x_tr,
                      theta=theta_tr,
                      f=f_tr,
                      method='PCGPwM')

pred_emu = emu_fresco.predict(x=x_test, theta=theta_test)

# Get the prediction means and variances
pred_m, pred_var = pred_emu.mean(), pred_emu.var()

# Check accuracy of emulator
SST = np.sum(np.square(f_test - np.mean(f_test.T, axis=1)))
SSE = np.sum(np.square(pred_m - f_test))
print('Rsq PCGP lin 1 = ', np.round(1 - SSE / SST, 2))

# Observed data
data = np.array([[26, 1243],
                 [31, 887.7],
                 [41, 355.5],
                 [51, 111.5],
                 [61, 26.5],
                 [71, 10.4],
                 [76, 8.3],
                 [81, 7.3],
                 [91, 17.2],
                 [101, 37.6],
                 [111, 48.7],
                 [121, 38.9],
                 [131, 32.4],
                 [141, 36.4],
                 [151, 61.9]])

x = data[:, 0].reshape((15, 1))
y = data[:, 1].reshape((15, 1))

obsvar = (0.2 * y) ** 2

# Fit calibrator with different methods and samplers
cal_1 = calibrator(emu=emu_fresco,
                   y=y,
                   x=x,
                   thetaprior=prior_fresco_Ca,
                   method='directbayes',
                   yvar=obsvar)

visualize(cal_1, emu_fresco, x_tr, theta_tr)

# cal_2 = calibrator(emu=emu_fresco,
#                     y=y,
#                     x=x,
#                     thetaprior=prior_fresco_Ca, 
#                     method='directbayeswoodbury',
#                     yvar=obsvar)

# visualize(cal_2, emu_fresco, x_tr, theta_tr)

# cal_3 = calibrator(emu=emu_fresco,
#                     y=y,
#                     x=x,
#                     thetaprior=prior_fresco_Ca, 
#                     method='directbayeswoodbury',
#                     yvar=obsvar,
#                     args={'sampler': 'PTLMC',
#                           'usedir': True})

# visualize(cal_3, emu_fresco, x_tr, theta_tr)

# plt.plot(emu_fresco.predict(x=x, theta=parameter_values).mean())
# plt.plot(y)


