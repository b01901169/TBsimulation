"""
==========================================================================
Illustration of prior and posterior Gaussian process for different kernels
==========================================================================

This example illustrates the prior and posterior of a GPR with different
kernels. Mean, standard deviation, and 10 samples are shown for both prior
and posterior.
"""
print(__doc__)

# Authors: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#
# License: BSD 3 clause

import numpy as np

from matplotlib import pyplot as plt

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import (RBF, Matern, RationalQuadratic,
                                              ExpSineSquared, DotProduct,
                                              ConstantKernel, WhiteKernel)

from kernel import composedKernel


if __name__ == "__main__":

    f1 = lambda X: np.array([1 for i in range(len(X))])
    f2 = lambda X: np.array([(float(X[i])+10)*0.2 for i in range(len(X))])
    f3 = lambda X: np.array([(float((X[i])-10)/3)**2 for i in range(len(X))])

    k1 = RBF(length_scale=0.6)
    k2 = RBF(length_scale=0.8)
    k3 = RBF(length_scale=1.0)
    k4 = RBF(length_scale=1.2)
    k5 = RBF(length_scale=1.4)
    k6 = RBF(length_scale=1.6)
    k7 = RBF(length_scale=1.8)
    k8 = RBF(length_scale=2.0)
    k9 = RBF(length_scale=2.2)
    k10 = RBF(length_scale=2.4)
    
    #X = np.array([[1,2,3,4]])
    #Y = np.array([[4,3,2,1]])

    # kernels = [1.0 * RBF(length_scale=1.0, length_scale_bounds=(1e-1, 10.0)),
    #            1.0 * RationalQuadratic(length_scale=1.0, alpha=0.1),
    #            1.0 * ExpSineSquared(length_scale=1.0, periodicity=3.0,
    #                                 length_scale_bounds=(0.1, 10.0),
    #                                 periodicity_bounds=(1.0, 10.0)),
    #            ConstantKernel(0.1, (0.01, 10.0))
    #                * (DotProduct(sigma_0=1.0, sigma_0_bounds=(0.0, 10.0)) ** 2),
    #            1.0 * Matern(length_scale=1.0, length_scale_bounds=(1e-1, 10.0),
    #                         nu=1.5)]
    
    time_horizon = 2
    fig_index = 0
    grid_size = 1000
    bias_sample_size = 10
    upper_bound = 1
    posterior_sample_size = 0
    sub_visible_region = 3
    visible_region = 5
    tick_size = 30
    plot_detail = True
    #individual_alpha = 0.01
    random_seed = np.random.randint(0,10000)
    random_seed = 6935
    gp_alpha = 0.000001

    print("random seed: {0}".format(random_seed))
    np.random.seed(random_seed)

    # ============================== composed kernels ===================================

    coefficients = [f1] * time_horizon
    #coefficients = [lambda X: np.array([((1.0 - float(X[i])/(5.0 * upper_bound)))**(t/3.0 + 1) for i in range(len(X))]) for t in range(time_horizon) ]
    #subkernels = [RBF(length_scale=(np.random.randint(30,150))*0.001) for t in range(time_horizon)]
    #subkernels = [RBF(length_scale=0.02 + np.random.rand() * 0.08) + 0.5 * Matern(length_scale=0.20 + np.random.rand() * 0.30, nu=1 + np.random.rand()) for i in range(time_horizon)]
    subkernels = [RBF(length_scale=0.1) + 0.7 * Matern(length_scale=0.1, nu=1.5),
                  RBF(length_scale=0.04) + 0.6 * Matern(length_scale=0.15, nu=1.3),
                 ]

    
    ck = composedKernel(coefficients, subkernels)

    # ===================================================================================
    
    X_values = np.linspace(0, upper_bound, grid_size)
    X_ = X_values[:, np.newaxis]

    random_indices = np.random.randint(0, grid_size, bias_sample_size)
    #random_indices = [30, 51, 103, 301, 509, 650, 699, 923, 1032, 1099, 1130, 1193, 1240, 1430, 1590, 1602, 1790, 1838, 1902, 1984]
    random_X = X_[random_indices]

    y_whole_target = np.zeros((grid_size, 1))
    y_whole_bias = np.zeros((bias_sample_size, 1))

    y_whole_prior_mean = np.zeros((grid_size))
    y_whole_prior_variance = np.zeros((grid_size))
    y_whole_mean = np.zeros((grid_size))
    y_whole_variance = np.zeros((grid_size))

    whole_log_likelihood = 0

    coefficient_square_upper_bound_list = [np.max(np.square(coefficients[t](X_))) for t in range(time_horizon)]
    coefficient_square_upper_bound_sum = np.sum(coefficient_square_upper_bound_list)
    print("coefficient square upper bound sum: {0}".format(coefficient_square_upper_bound_sum))

    y_sample_list = []
    y_whole_posterior_samples = np.zeros((grid_size, posterior_sample_size))
    for i in range(time_horizon):
        gp = GaussianProcessRegressor(kernel=subkernels[i], optimizer=None, alpha=gp_alpha)

        # --------------------------- prior ------------------------------
        y_prior_mean, y_prior_std = gp.predict(X_, return_std=True)
    
        y_sample = gp.sample_y(X_, 1, np.random.randint(10000))
        #y_sample = gp.sample_y(X_, 1)
        y_sample_list.append(y_sample)
        y = y_sample[random_indices]

        # ------------------------- plot prior ---------------------------
        if plot_detail:
            plt.figure(i, figsize=(8, 8))
            #plt.subplot(2, 1, 1)
    
            #plt.plot(X_values, y_prior_mean, 'k', lw=3, zorder=9)
            #plt.fill_between(X_values, y_prior_mean - y_prior_std, y_prior_mean + y_prior_std,
            #             alpha=0.2, color='k')
    
            plt.plot(X_values, y_sample, lw=3)
            plt.scatter(random_X[:, 0], y, c='r', s=100, zorder=posterior_sample_size, edgecolors=(0, 0, 0))
            plt.xlim(0, upper_bound)
            #plt.ylim(-sub_visible_region, sub_visible_region)
            plt.ylim(-5, 3)
            #plt.title("Target: subfunction {0}".format(i+1), fontsize=12)
            plt.tick_params(axis="both", labelsize=tick_size, bottom=False, left=False, labelbottom=False, labelleft=False)
            plt.tight_layout()
            plt.savefig("figure/{0}_1_over_{1}.png".format(i, time_horizon))


        # ------------------------- posterior -----------------------------
        gp.fit(random_X, y)

        y_mean, y_std = gp.predict(X_, return_std=True)
        y_mean = y_mean[:,0]

        y_posterior_samples = gp.sample_y(X_, posterior_sample_size, np.random.randint(10000))
        #y_posterior_samples = y_posterior_samples[:,0,:]
        if posterior_sample_size:
            y_posterior_samples = y_posterior_samples[:,0,np.random.choice(posterior_sample_size, posterior_sample_size, replace=False)]
        else:
            y_posterior_samples = y_posterior_samples[:,0,:]

        # ---------------------- plot posterior ---------------------------
        if plot_detail:
            plt.figure(time_horizon+i, figsize=(8, 8))
            #plt.subplot(2, 1, 2)
    
            plt.plot(X_values, y_sample, lw=3)
            #plt.plot(X_values, y_mean, 'k', lw=3, zorder=9)
            plt.fill_between(X_values, y_mean - y_std, y_mean + y_std,
                             alpha=0.2, color='k')
    
            if posterior_sample_size:
                plt.plot(X_values, y_posterior_samples, lw=3)
            plt.scatter(random_X[:, 0], y, c='r', s=100, zorder=posterior_sample_size, edgecolors=(0, 0, 0))
            plt.xlim(0, upper_bound)
            #plt.ylim(-sub_visible_region, sub_visible_region)
            plt.ylim(-5, 3)
            #plt.title("Posterior: subfunction {0}".format(i+1), fontsize=12)
            plt.tick_params(axis="both", labelsize=tick_size, bottom=False, left=False, labelbottom=False, labelleft=False)
            plt.tight_layout()
            plt.savefig("figure/{0}_2_over_{1}.png".format(i, time_horizon))

        # -------------------- composition of subkernels ------------------------
        multiply_coefficient = coefficients[i](X_)
        y_whole_target += (y_sample[:,0] * multiply_coefficient)[:, np.newaxis]
        y_whole_posterior_samples += (np.repeat(multiply_coefficient[:,np.newaxis], posterior_sample_size, axis=1) * y_posterior_samples)

        y_whole_prior_mean += y_prior_mean * multiply_coefficient
        y_whole_prior_variance += (y_prior_std * multiply_coefficient)**2
        y_whole_mean += y_mean * multiply_coefficient
        y_whole_variance += (y_std * multiply_coefficient)**2

        y_whole_bias += (y * multiply_coefficient[random_indices, np.newaxis])
        whole_log_likelihood += gp.log_marginal_likelihood(gp.kernel_.theta)

    y_whole_prior_std = np.sqrt(y_whole_prior_variance)
    y_whole_std = np.sqrt(y_whole_variance)

    # ================= decomposed Gaussian process regression ==================
    plt.figure(time_horizon*2, figsize=(24,12))
    #plt.subplot(2, 1, 1)

    #plt.plot(X_values, y_whole_prior_mean, 'k', lw=3, zorder=9)
    #plt.fill_between(X_values, y_whole_prior_mean - y_whole_prior_std, y_whole_prior_mean + y_whole_prior_std,
    #                 alpha=0.2, color='k')

    plt.plot(X_values, y_whole_target, lw=3)
    plt.scatter(random_X[:, 0], y_whole_bias, c='r', s=100, zorder=posterior_sample_size, edgecolors=(0, 0, 0))
    plt.xlim(0, upper_bound)
    #plt.ylim(-visible_region, visible_region)
    plt.ylim(-6, 3)
    #plt.title("Target", fontsize=12)
    plt.tick_params(axis="both", labelsize=tick_size, bottom=False, left=False, labelbottom=False, labelleft=False)
    plt.tight_layout()
    plt.savefig("figure/entire_prior.png")

    plt.figure(time_horizon*2+1, figsize=(24,12))
    #plt.subplot(2, 1, 2)
    plt.plot(X_values, y_whole_target, lw=3)
    #plt.plot(X_values, y_whole_mean, 'k', lw=3, zorder=9)
    plt.fill_between(X_values, y_whole_mean - y_whole_std, y_whole_mean + y_whole_std,
                     alpha=0.2, color='g')

    if posterior_sample_size:
        plt.plot(X_values, y_whole_posterior_samples, lw=3)
    plt.scatter(random_X[:, 0], y_whole_bias, c='r', s=100, zorder=posterior_sample_size, edgecolors=(0, 0, 0))
    plt.xlim(0, upper_bound)
    #plt.ylim(-visible_region, visible_region)
    plt.ylim(-6, 3)
    #plt.title("Posterior (decomposed GP regression)\n avearage std: %.3f"
    #          % (np.mean(y_whole_std)),
    #          fontsize=12)
    plt.tick_params(axis="both", labelsize=tick_size, bottom=False, left=False, labelbottom=False, labelleft=False)
    plt.tight_layout()
    plt.savefig("figure/decomposedGPs.png")


    # =========== without decomposition Gaussian process regression =============
    gp = GaussianProcessRegressor(kernel=ck, optimizer=None, alpha=gp_alpha*coefficient_square_upper_bound_sum)

    y_mean, y_std = gp.predict(X_, return_std=True)
    y_samples = gp.sample_y(X_, posterior_sample_size, np.random.randint(10000))

    # ------ Plot prior ------
    plt.figure(time_horizon*2+2, figsize=(12,12))

    #plt.subplot(2, 1, 1)

    plt.plot(X_values, y_mean, 'k', lw=3, zorder=9)
    plt.fill_between(X_values, y_mean - y_std, y_mean + y_std,
                     alpha=0.2, color='k')
    if posterior_sample_size:
        plt.plot(X_values, y_samples, lw=1)
    plt.xlim(0, upper_bound)
    #plt.ylim(-visible_region, visible_region)
    plt.ylim(-6, 3)
    #plt.title("Prior (kernel:  %s)" % gp.kernel, fontsize=12)
    plt.tick_params(axis="both", labelsize=tick_size, bottom=False, left=False, labelbottom=False, labelleft=False)
    plt.tight_layout()

    # -------------------------- Generate data and fit GP ------------------------
    gp.fit(random_X, y_whole_bias)

    y_mean, y_std = gp.predict(X_, return_std=True)
    y_mean = y_mean[:,0]

    y_samples = gp.sample_y(X_, posterior_sample_size, np.random.randint(10000))
    y_samples = y_samples[:,0,:]

    # ------ Plot posterior ------
    plt.figure(time_horizon*2+3, figsize=(24,12))
    #plt.subplot(2, 1, 2)

    plt.plot(X_values, y_whole_target, lw=3)
    #plt.plot(X_values, y_mean, 'k', lw=3, zorder=9)
    plt.fill_between(X_values, y_mean - y_std, y_mean + y_std,
                     alpha=0.2, color='b')

    if posterior_sample_size:
        plt.plot(X_values, y_samples, lw=1)
    plt.scatter(random_X[:, 0], y_whole_bias, c='r', s=100, zorder=posterior_sample_size, edgecolors=(0, 0, 0))
    plt.xlim(0, upper_bound)
    #plt.ylim(-visible_region, visible_region)
    plt.ylim(-6, 3)
    #plt.title("Posterior (GP regression)\n average std: %.3f" % (np.mean(y_std)), fontsize=12)
    plt.tick_params(axis="both", labelsize=tick_size, bottom=False, left=False, labelbottom=False, labelleft=False)
    plt.tight_layout()
    plt.savefig("figure/entireGP.png")
    
    # ------ Plot posterior comparison ------
    plt.figure(time_horizon*2+4, figsize=(24,12))
    #plt.subplot(2, 1, 2)

    plt.plot(X_values, y_whole_target, lw=3)
    #plt.plot(X_values, y_mean, 'k', lw=3, zorder=9)
    plt.fill_between(X_values, y_whole_mean - y_whole_std, y_whole_mean + y_whole_std,
                     alpha=0.2, color='g', label="Decomposed Gaussian Process Regression")
    plt.fill_between(X_values, y_mean - y_std, y_mean + y_std,
                     alpha=0.2, color='b', label="Gaussian Process Regression")

    if posterior_sample_size:
        plt.plot(X_values, y_samples, lw=1)
    plt.scatter(random_X[:, 0], y_whole_bias, c='r', s=100, zorder=posterior_sample_size, edgecolors=(0, 0, 0))
    plt.xlim(0, upper_bound)
    #plt.ylim(-visible_region, visible_region)
    plt.ylim(-6, 3)
    #plt.title("Posterior (GP regression)\n Average std: %.3f, %.3f" % (np.mean(y_std), np.mean(y_whole_std)), fontsize=12)
    print("Posterior (GP regression)\n Average std: %.3f, %.3f" % (np.mean(y_std), np.mean(y_whole_std)))
    plt.legend(fontsize=50, loc="lower right")
    plt.tick_params(axis="both", labelsize=tick_size, bottom=False, left=False, labelbottom=False, labelleft=False)
    plt.tight_layout()
    plt.savefig("figure/comparison.png")

    #plt.show()

