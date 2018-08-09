import numpy as np
import scipy
from sklearn.gaussian_process.kernels import (RBF, Matern, RationalQuadratic,
                                              ExpSineSquared, DotProduct,
                                              ConstantKernel, WhiteKernel)

import pandas

from decomposition import *

if __name__ == "__main__":
    date = "0809"
    J = 5
    # a = 1
    # b = 1
    upper_bound = 1
    grid_size = 1000
    gp_alpha = 0.00001
    dimension = 1
    random_seed = np.random.randint(1000)
    constraints = None
    discrete = False

    X_ = np.reshape(np.linspace(0, upper_bound, grid_size), (grid_size, 1))


    #kernelList = [RBF(length_scale=0.2), RBF(length_scale=0.5)]
    #kernelList = [Matern(length_scale=0.1*i) for i in range(1,J+1)]
    kernelList = [RBF(length_scale=0.05*i) for i in range(1,J+1)]

    kernel = sum(kernelList)
    g = lambda x: np.sum(x)

    targetList = [GaussianProcessRegressor(kernel=kernelList[i], optimizer=None, alpha=gp_alpha).sample_y(X_, 1, random_seed) for i in range(J)]
    fList = [randomize(smoothify(targetList[i], X_), gp_alpha) for i in range(J)]

    function_bounds = np.zeros(J)
    for i in range(J):
        function_bounds[i] = np.mean(np.prod([np.abs(targetList[j]) for j in np.delete(np.arange(J), i)], axis=0))

    b_scale = np.sum(function_bounds)
    gList = [lambda x: function_bounds[i] for i in range(J)]

    decomposition = Decomposition(J, fList, g, gList, kernelList)

    max_derivative_list = [maxDerivative(targetList[i], grid_size) for i in range(J)]

    # -------------------------------- experiment -----------------------------------
    total_count = 1
    total_run = 200
    a_count = 8
    a_list = np.array([0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.5, 1]) * np.max(max_derivative_list)
    b_count = 8 
    b_list = np.array([0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.5, 1])

    GPUCB_scores = np.zeros((a_count, b_count, total_count))
    decomposedGPUCB_scores = np.zeros((a_count, b_count, total_count))
    for a_index in range(a_count):
        a = a_list[a_index]
        for b_index in range(b_count):
            b = b_list[b_index]
            for count in range(total_count):
                GPUCBsolver = GPUCB(decomposition.get_function_value, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, a=a, b=b, X_=X_, discrete=discrete, linear=True) # linear arg only changes the beta_t used in exploration
                GPUCBsolver.run(total_run)
                GPUCB_scores[a_index, b_index, count] = GPUCBsolver.regret

                decomposedGPUCBsolver = DecomposedGPUCB(decomposition, kernelList, dimension, upper_bound, constraints, gp_alpha=gp_alpha, a=a, b=b, X_=X_, discrete=discrete)
                decomposedGPUCBsolver.run(total_run)
                decomposedGPUCB_scores[a_index, b_index, count] = decomposedGPUCBsolver.regret

    GPUCB_df = pandas.DataFrame(data=GPUCB_scores[:,:,0], columns=b_list, index=a_list)
    decomposedGPUCB_df = pandas.DataFrame(data=decomposedGPUCB_scores[:,:,0], columns=b_list, index=a_list)

    GPUCB_df.to_csv(path_or_buf='result/GPUCB_result_{0}.csv'.format(date))
    decomposedGPUCB_df.to_csv(path_or_buf='result/decomposedGPUCB_result_{0}.csv'.format(date))

