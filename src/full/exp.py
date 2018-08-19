import numpy as np
import scipy
from sklearn.gaussian_process.kernels import (RBF, Matern, RationalQuadratic,
                                              ExpSineSquared, DotProduct,
                                              ConstantKernel, WhiteKernel)

import pandas as pd
import argparse
import pickle

from decomposition import *

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Decomposed GPUCB and GPUCB comparison')
    parser.add_argument('-n', '--name', required=True, help='Input the name to save')

    output_path = "synthetic/linear/"

    args = parser.parse_args()

    filename = args.name
    J = 5
    # a = 1
    # b = 1
    upper_bound = 1
    grid_size = 1000
    gp_alpha_list = [0.01] * J
    gp_alpha = 0.01 * J
    delta = 0.05
    dimension = 1
    random_seed = np.random.randint(1000)
    constraints = None
    discrete = True

    X_ = np.reshape(np.linspace(0, upper_bound, grid_size), (grid_size, 1))


    #kernelList = [RBF(length_scale=0.2), RBF(length_scale=0.5)]
    #kernelList = [Matern(length_scale=0.1*i) for i in range(1,J+1)]
    kernelList = [RBF(length_scale=0.005*i) for i in range(1,J+1)]

    kernel = sum(kernelList)
    g = lambda x: np.sum(x)

    targetList = [GaussianProcessRegressor(kernel=kernelList[i], optimizer=None, alpha=gp_alpha).sample_y(X_, 1, random_seed) for i in range(J)]
    fList = [randomizify(smoothify(targetList[i], X_), gp_alpha) for i in range(J)]

    function_bounds = np.zeros(J)
    for i in range(J):
        function_bounds[i] = np.mean(np.prod([np.abs(targetList[j]) for j in np.delete(np.arange(J), i)], axis=0))

    b_scale = np.sum(function_bounds)
    # gList = [lambda x: function_bounds[i] for i in range(J)]
    gList = [lambda x: 1.0 for i in range(J)]

    decomposition = Decomposition(J, fList, g, gList, kernelList)

    max_derivative_list = [maxDerivative(targetList[i], grid_size) for i in range(J)]

    # -------------------------------- experiment -----------------------------------
    total_count = 10
    total_run = 300
    a_count = 5
    #a_list = np.array([0.002]) * np.mean(max_derivative_list)
    a_list = np.array([0.005, 0.01, 0.02, 0.05, 0.1]) * np.mean(max_derivative_list)
    b_count = 5
    #b_list = np.array(np.arange(0.05, 0.06, 0.01))
    b_list = np.array(np.arange(0.01, 0.06, 0.01))

    GPUCB_scores = np.zeros((a_count, b_count))
    decomposedGPUCB_scores = np.zeros((a_count, b_count))
    GPUCB_regret_list = np.zeros((a_count, b_count, total_run))
    decomposed_regret_list = np.zeros((a_count, b_count, total_run))

    for count in range(total_count):
        # GPUCB_scores = np.zeros((a_count, b_count))
        # decomposedGPUCB_scores = np.zeros((a_count, b_count))

        for a_index in range(a_count):
            a = a_list[a_index]
            for b_index in range(b_count):
                b = b_list[b_index]

                initial_point = X_[np.random.randint(grid_size)]

                print ("\nGPUCB count:{0}, a index:{1}, b index:{2}...".format(count, a_index, b_index))
                GPUCBsolver = GPUCB(decomposition.get_function_value, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, delta=delta, a=a, b=b, X_=X_, initial_point=initial_point, discrete=discrete, linear=True) # linear arg only changes the beta_t used in exploration
                GPUCBsolver.run(total_run)
                GPUCB_scores[a_index, b_index] += GPUCBsolver.regret
                GPUCB_regret_list[a_index, b_index] += np.array(GPUCBsolver.regret_list)

                print ("\ndecomposed count:{0}, a index:{1}, b index:{2}...".format(count, a_index, b_index))
                decomposedGPUCBsolver = DecomposedGPUCB(decomposition, kernelList, dimension, upper_bound, constraints, gp_alpha=gp_alpha_list, delta=delta, a=a, b=b, X_=X_, initial_point=initial_point, discrete=discrete)
                decomposedGPUCBsolver.run(total_run)
                decomposedGPUCB_scores[a_index, b_index] += decomposedGPUCBsolver.regret
                decomposed_regret_list[a_index, b_index] += np.array(decomposedGPUCBsolver.regret_list)

    GPUCB_df = pd.DataFrame(data=GPUCB_scores, columns=b_list, index=a_list)
    decomposedGPUCB_df = pd.DataFrame(data=decomposedGPUCB_scores, columns=b_list, index=a_list)

    GPUCB_df.to_csv(path_or_buf=output_path+'GPUCB_result_{0}.csv'.format(filename))
    decomposedGPUCB_df.to_csv(path_or_buf=output_path+'decomposedGPUCB_result_{0}.csv'.format(filename))

    decomposed_regret_list /= total_count
    GPUCB_regret_list /= total_count
    pickle.dump((GPUCB_regret_list, decomposed_regret_list), open(output_path+"regret_list_{0}.p".format(filename), 'wb'))

