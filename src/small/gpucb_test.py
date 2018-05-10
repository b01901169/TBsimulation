"""
==========================================================================
Illustration of prior and posterior Gaussian process for different kernels
==========================================================================

This example illustrates the prior and posterior of a GPR with different
kernels. Mean, standard deviation, and 10 samples are shown for both prior
and posterior.
"""
print(__doc__)

import numpy as np

from matplotlib import pyplot as plt

import scipy
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import (RBF, Matern, RationalQuadratic,
                                              ExpSineSquared, DotProduct,
                                              ConstantKernel, WhiteKernel)

from kernel import composedKernel


def decomposed_gpucb(coefficients, subkernels, target_function_list, X_, sample_X_indices, iterations, gp_alpha, upper_bound, f_output=None):
    assert(len(coefficients) == len(subkernels))
    dimension = X_.shape[1]
    time_horizon = len(coefficients)
    a = 1.0
    b = 1.0
    delta = 0.01
    sample_size, feature_size = X_.shape
    current_sample_size = len(sample_X_indices)

    target_function = np.zeros(sample_size)
    for t in range(time_horizon):
        sub_target_function = target_function_list[t]
        multiply_coefficient = coefficients[t](X_)
        target_function += sub_target_function[:,0] * multiply_coefficient

    best_x = np.argmax(target_function)
    maximum_value = target_function[best_x]
    minimum_value = np.min(target_function)
    print("Global maximum: {0}, global minimum: {1}, range: {2} ".format(maximum_value, minimum_value, maximum_value - minimum_value))

    if f_output:
        f_output.write("decomposed gpucb, global max, {0}, global min, {1}, range, {2}, ".format(maximum_value, minimum_value, maximum_value - minimum_value))

    existing_regret = np.mean(maximum_value - target_function[sample_X_indices])
    average_regret_list = [existing_regret]
    for iteration in range(current_sample_size, iterations):
        #gp_list = []
        sample_X = X_[sample_X_indices]
        y_mean = np.zeros(sample_size)
        y_variance = np.zeros(sample_size)
        for t in range(time_horizon):
            sub_target_function = target_function_list[t]
            sample_y = sub_target_function[sample_X_indices]
            gp = GaussianProcessRegressor(kernel=subkernels[t], optimizer=None, alpha=gp_alpha)
            #gp = GaussianProcessRegressor(kernel=subkernels[t], alpha=gp_alpha)
            gp.fit(sample_X, sample_y)
            y_sub_mean, y_sub_std = gp.predict(X_, return_std=True)

            #multiply_coefficient = np.sqrt(np.abs(coefficients[t](X_)))
            multiply_coefficient = coefficients[t](X_)
            #multiply_coefficient = np.square(coefficients[t](X_))
            #multiply_coefficient = 1

            #print(y_sub_mean.shape)
            #print(y_sub_std.shape)
            y_mean += y_sub_mean[:,0] * multiply_coefficient
            y_variance += (y_sub_std * multiply_coefficient)**2

        y_std = np.sqrt(y_variance)

        b_iteration = 2 * np.log(2 * iteration**2 * np.pi**2 * time_horizon / (3 * delta)) + 2 * dimension * np.log(iteration**2 * dimension * b * upper_bound * np.sqrt(np.log(4 * dimension * a / delta)))
        #b_iteration = 2 * np.log(2 * iteration**2 * np.pi**2 / (3 * delta)) + 2 * dimension * np.log(iteration**2 * dimension * b * upper_bound * np.sqrt(np.log(4 * dimension * a / delta)))

        ucb_y = y_mean + np.sqrt(b_iteration) * y_std
        argmax_x = np.argmax(ucb_y)

        sample_X_indices = np.concatenate((sample_X_indices, [argmax_x]), axis=0)

        average_regret = np.mean(maximum_value - target_function[sample_X_indices])
        average_regret_list.append(average_regret)
        current_maximum = np.max(target_function[sample_X_indices])
        print("Decomposed GP, iteration: {0}, average regret: {1}, current maximum: {2}, b_t: {3}".format(iteration, average_regret, current_maximum, b_iteration))

    if f_output:
        f_output.write(", ".join([str(x) for x in average_regret_list]) + "\n")

    return sample_X_indices

def entire_gpucb(coefficients, subkernels, target_function_list, X_, sample_X_indices, iterations, gp_alpha, upper_bound, f_output=None):
    assert(len(coefficients) == len(subkernels))
    dimension = X_.shape[1]
    time_horizon = len(coefficients)
    a = 1.0 
    b = 1.0 * time_horizon
    delta = 0.01
    sample_size, feature_size = X_.shape
    current_sample_size = len(sample_X_indices)

    target_function = np.zeros(sample_size)
    for t in range(time_horizon):
        sub_target_function = target_function_list[t]
        multiply_coefficient = coefficients[t](X_)
        target_function += sub_target_function[:,0] * multiply_coefficient

    best_x = np.argmax(target_function)
    maximum_value = target_function[best_x]
    minimum_value = np.min(target_function)
    print("Global maximum: {0}, global minimum: {1}, range: {2}".format(maximum_value, minimum_value, maximum_value - minimum_value))

    if f_output:
        f_output.write("entire gpucb, global max, {0}, global min, {1}, range, {2}, ".format(maximum_value, minimum_value, maximum_value - minimum_value))

    existing_regret = np.mean(maximum_value - target_function[sample_X_indices])
    average_regret_list = [existing_regret]

    # ------------------ composition of subkernels --------------------
    ck = composedKernel(coefficients, subkernels)
    #ck = RBF(alpha=0.001)

    # ---------------------------- GPUCB ------------------------------
    for iteration in range(current_sample_size, iterations):
        #gp_list = []
        sample_X = X_[sample_X_indices]
        sample_y = target_function[sample_X_indices]

        # ----------------------- gp regression -----------------------
        gp = GaussianProcessRegressor(kernel=ck, optimizer=None, alpha=gp_alpha)
        gp.fit(sample_X, sample_y)

        y_mean, y_std = gp.predict(X_, return_std=True)

        # ---------------------- arg-maximization ---------------------
        b_iteration = 2 * np.log(2 * iteration**2 * np.pi**2 / (3 * delta)) + 2 * dimension * np.log(iteration**2 * dimension * b * upper_bound * np.sqrt(np.log(4 * dimension * a / delta)) )
        ucb_y = y_mean + np.sqrt(b_iteration) * y_std
        argmax_x = np.argmax(ucb_y)

        sample_X_indices = np.concatenate((sample_X_indices, [argmax_x]), axis=0)

        average_regret = np.mean(maximum_value - target_function[sample_X_indices])        
        average_regret_list.append(average_regret)
        current_maximum = np.max(target_function[sample_X_indices])
        print("Entire GP, iteration: {0}, average regret: {1}, current maximum: {2}, b_t: {3}".format(iteration, average_regret, current_maximum, b_iteration))

    if f_output:
        f_output.write(", ".join([str(x) for x in average_regret_list]) + "\n")

    return sample_X_indices

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
    
    time_horizon = 32
    fig_index = 0
    grid_size = 1000
    bias_sample_size = 1
    upper_bound = 1
    posterior_sample_size = 30
    visible_region = 5
    plot_detail = False
    individual_alpha = 0.01
    gp_alpha = 1e-10
    iterations = 100
    output_filename = "result/result_0508.csv"

    f_output = open(output_filename, "a")
    random_seed = np.random.randint(0,10000)

    np.random.seed()

    # ======================== coefficients and sub-kernels ===========================
    coefficients = [lambda X: np.array([((1 - float(X[i])/(5.0 * upper_bound)))**(t/3.0 + 1) for i in range(len(X))]) for t in range(time_horizon) ]
    #coefficients = [f1] * time_horizon
    subkernels = [RBF(length_scale=(np.random.randint(20,50))*0.001,) + WhiteKernel(noise_level=individual_alpha) for i in range(time_horizon)]
    #subkernels = [RBF(length_scale=(5)*0.005,) + WhiteKernel(noise_level=individual_alpha) for i in range(time_horizon)]
    #subkernels = [Matern(length_scale=(i+10)*0.05, nu=1.5) + WhiteKernel(noise_level=individual_alpha) for i in range(time_horizon)]
    
    ck = composedKernel(coefficients, subkernels)

    # =================================================================================
 
    X_values = np.linspace(0, upper_bound, grid_size)
    X_ = X_values[:, np.newaxis]

    random_indices = np.random.randint(0, grid_size, bias_sample_size)
    random_X = X_[random_indices]

    y_whole_target = np.zeros((grid_size, 1))
    y_whole_bias = np.zeros((bias_sample_size, 1))

    y_whole_prior_mean = np.zeros((grid_size))
    y_whole_prior_variance = np.zeros((grid_size))
    y_whole_mean = np.zeros((grid_size))
    y_whole_variance = np.zeros((grid_size))

    whole_log_likelihood = 0

    y_sample_list = []
    y_whole_posterior_samples = np.zeros((grid_size, posterior_sample_size))
    for i in range(time_horizon):
        gp = GaussianProcessRegressor(kernel=subkernels[i], optimizer=None, alpha=gp_alpha)

        # --------------------------- prior ------------------------------
        y_prior_mean, y_prior_std = gp.predict(X_, return_std=True)
    
        y_sample = gp.sample_y(X_, 1, random_seed)
        y_sample_list.append(y_sample)

        """
        # ------------------------- plot prior ---------------------------
        if plot_detail:
            plt.figure(i, figsize=(8, 8))
            plt.subplot(2, 1, 1)
    
            plt.plot(X_values, y_prior_mean, 'k', lw=3, zorder=9)
            plt.fill_between(X_values, y_prior_mean - y_prior_std, y_prior_mean + y_prior_std,
                         alpha=0.2, color='k')
    
            plt.plot(X_values, y_sample, lw=1)
            plt.xlim(0, upper_bound)
            plt.ylim(-3, 3)
            plt.title("Target (kernel:  %s)" % subkernels[i], fontsize=12)

        # ------------------------- posterior -----------------------------
        y = y_sample[random_indices]
        gp.fit(random_X, y)

        y_mean, y_std = gp.predict(X_, return_std=True)
        y_mean = y_mean[:,0]

        y_posterior_samples = gp.sample_y(X_, posterior_sample_size)
        #y_posterior_samples = y_posterior_samples[:,0,:]
        y_posterior_samples = y_posterior_samples[:,0,np.random.choice(posterior_sample_size, posterior_sample_size, replace=False)]

        # ---------------------- plot posterior ---------------------------
        if plot_detail:
            plt.subplot(2, 1, 2)
    
            plt.plot(X_values, y_mean, 'k', lw=3, zorder=9)
            plt.fill_between(X_values, y_mean - y_std, y_mean + y_std,
                             alpha=0.2, color='k')
    
            plt.plot(X_values, y_posterior_samples, lw=1)
            plt.scatter(random_X[:, 0], y, c='r', s=50, zorder=posterior_sample_size, edgecolors=(0, 0, 0))
            plt.xlim(0, upper_bound)
            plt.ylim(-3, 3)
            plt.title("Posterior (kernel: %s)\n Log-Likelihood: %.3f"
                      % (gp.kernel_, gp.log_marginal_likelihood(gp.kernel_.theta)),
                      fontsize=12)
            plt.tight_layout()
        #"""

        # -------------------- composition of subkernels ------------------------
        multiply_coefficient = coefficients[i](X_)
        y_whole_target += (y_sample[:,0] * multiply_coefficient)[:, np.newaxis]

        y_whole_prior_mean += y_prior_mean * multiply_coefficient
        y_whole_prior_variance += (y_prior_std * multiply_coefficient)**2

        #y_whole_posterior_samples += (np.repeat(multiply_coefficient[:,np.newaxis], posterior_sample_size, axis=1) * y_posterior_samples)
        #y_whole_mean += y_mean * multiply_coefficient
        #y_whole_variance += (y_std * multiply_coefficient)**2

        #y_whole_bias += (y * multiply_coefficient[random_indices, np.newaxis])
        #whole_log_likelihood += gp.log_marginal_likelihood(gp.kernel_.theta)

    y_whole_prior_std = np.sqrt(y_whole_prior_variance)
    y_whole_std = np.sqrt(y_whole_variance)

    decomposed_gpucb_indices = decomposed_gpucb(coefficients, subkernels, y_sample_list, X_, random_indices, iterations, gp_alpha, upper_bound, f_output)
    entire_gpucb_indices = entire_gpucb(coefficients, subkernels, y_sample_list, X_, random_indices, iterations, gp_alpha, upper_bound, f_output)

    f_output.close()

    """
    # ================= decomposed Gaussian process regression ==================
    plt.figure(time_horizon, figsize=(12,12))
    plt.subplot(2, 1, 1)

    plt.plot(X_values, y_whole_prior_mean, 'k', lw=3, zorder=9)
    plt.fill_between(X_values, y_whole_prior_mean - y_whole_prior_std, y_whole_prior_mean + y_whole_prior_std,
                     alpha=0.2, color='k')

    plt.plot(X_values, y_whole_target, lw=1)
    plt.xlim(0, upper_bound)
    plt.ylim(-visible_region, visible_region)
    plt.title("Target (composition of subkernels)", fontsize=12)

    plt.subplot(2, 1, 2)
    plt.plot(X_values, y_whole_mean, 'k', lw=3, zorder=9)
    plt.fill_between(X_values, y_whole_mean - y_whole_std, y_whole_mean + y_whole_std,
                     alpha=0.2, color='k')

    plt.plot(X_values, y_whole_posterior_samples, lw=1)
    plt.scatter(random_X[:, 0], y_whole_bias, c='r', s=50, zorder=posterior_sample_size, edgecolors=(0, 0, 0))
    plt.xlim(0, upper_bound)
    plt.ylim(-visible_region, visible_region)
    plt.title("Posterior (kernel: composition of subkernels)\n avearage std: %.3f"
              % (np.mean(y_whole_std)),
              fontsize=12)
    plt.tight_layout()


    # =========== without decomposition Gaussian process regression =============
    gp = GaussianProcessRegressor(kernel=ck, optimizer=None, alpha=individual_alpha*time_horizon)

    y_mean, y_std = gp.predict(X_, return_std=True)
    y_samples = gp.sample_y(X_, posterior_sample_size)

    # Plot prior
    plt.figure(time_horizon+1, figsize=(12,12))

    plt.subplot(2, 1, 1)

    plt.plot(X_values, y_mean, 'k', lw=3, zorder=9)
    plt.fill_between(X_values, y_mean - y_std, y_mean + y_std,
                     alpha=0.2, color='k')
    plt.plot(X_values, y_samples, lw=1)
    plt.xlim(0, upper_bound)
    plt.ylim(-visible_region, visible_region)
    plt.title("Prior (kernel:  %s)" % ck, fontsize=12)

    # -------------------------- Generate data and fit GP ------------------------
    gp.fit(random_X, y_whole_bias)

    y_mean, y_std = gp.predict(X_, return_std=True)
    y_mean = y_mean[:,0]

    y_samples = gp.sample_y(X_, posterior_sample_size)
    y_samples = y_samples[:,0,:]

    # Plot posterior
    plt.subplot(2, 1, 2)

    plt.plot(X_values, y_mean, 'k', lw=3, zorder=9)
    plt.fill_between(X_values, y_mean - y_std, y_mean + y_std,
                     alpha=0.2, color='k')

    plt.plot(X_values, y_samples, lw=1)
    plt.scatter(random_X[:, 0], y_whole_bias, c='r', s=50, zorder=posterior_sample_size, edgecolors=(0, 0, 0))
    plt.xlim(0, upper_bound)
    plt.ylim(-visible_region, visible_region)
    plt.title("Posterior (kernel: %s)\n average std: %.3f"
              % (gp.kernel_, np.mean(y_std)),
              fontsize=12)
    plt.tight_layout()
    
    plt.show()
    #"""
