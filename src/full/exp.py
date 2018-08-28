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
    parser.add_argument('-s', '--scale_down', required=True, help='Input the scale down factor')
    parser.add_argument('-iteration', '--iteration', default=300, help='Input the total iterations')
    parser.add_argument('-count', '--count', default=10, help='Input the total count')
    parser.add_argument('-a', '--a', required=True, help='Input the a value')
    parser.add_argument('-b', '--b', required=True, help='Input the b value')
    parser.add_argument('-J', '--J', required=True, help='Input the J value')

    args = parser.parse_args()
    scale_down_factor = float(args.scale_down)
    filename = "{0}_scale{1}_a{2}_b{3}".format(args.name, args.scale_down, args.a, args.b)
    total_run = int(args.iteration)
    total_count = int(args.count)

    output_path = "synthetic/linear/"

    J = int(args.J)
    # a = 1
    # b = 1
    upper_bound = 1
    grid_size = 1000
    gp_alpha_list = [0.01] * J
    gp_alpha = 0.01 * J
    delta = 0.05
    dimension = 1
    constraints = None
    discrete = True
    maxmin = "max"

    X_ = np.reshape(np.linspace(0, upper_bound, grid_size), (grid_size, 1))


    #kernelList = [RBF(length_scale=0.2), RBF(length_scale=0.5)]
    #kernelList = [Matern(length_scale=0.15*i) for i in range(1,J+1)]
    kernelList = [RBF(length_scale=0.005*i) for i in range(1,J+1)]

    kernel = sum(kernelList)
    g = lambda x: np.sum(x)

    # ========================== experimental design ========================
    # ================================ experimental records ===============================
    GPUCB_scores = np.zeros(total_count)
    decomposedGPUCB_scores = np.zeros(total_count)
    EI_scores = np.zeros(total_count)
    POI_scores = np.zeros(total_count)
    decomposedEI_scores = np.zeros(total_count)
    decomposedPOI_scores = np.zeros(total_count)

    GPUCB_regret_list = np.zeros((total_count, total_run))
    decomposedGPUCB_regret_list = np.zeros((total_count, total_run))
    EI_regret_list = np.zeros((total_count, total_run))
    POI_regret_list = np.zeros((total_count, total_run))
    decomposedEI_regret_list = np.zeros((total_count, total_run))
    decomposedPOI_regret_list = np.zeros((total_count, total_run))

    f_output = open(output_path + "scores_report_{0}.csv".format(filename), 'w')
    f_output.write("round, GPUCB, decomposed GPUCB, EI, decomposed EI, POI, decomposed POI")
    for count in range(total_count):
        # ============================ generating functions ========================
        random_seed = np.random.randint(10000)
        print("random seed: {0}".format(random_seed))
        targetList = [GaussianProcessRegressor(kernel=kernelList[i], optimizer=None, alpha=gp_alpha).sample_y(X_, 1, random_state=random_seed) for i in range(J)]
        fList = [randomizify(smoothify(targetList[i], X_), gp_alpha) for i in range(J)]

        function_bounds = np.zeros(J)
        for i in range(J):
            function_bounds[i] = np.mean(np.prod([np.abs(targetList[j]) for j in np.delete(np.arange(J), i)], axis=0))

        b_scale = np.sum(function_bounds)
        # gList = [lambda x: function_bounds[i] for i in range(J)]
        gList = [lambda x: 1.0 for i in range(J)]

        decomposition = Decomposition(J, fList, g, gList)
        max_derivative_list = [maxDerivative(targetList[i], grid_size) for i in range(J)]

        a = float(args.a)
        b = float(args.b) * np.max(max_derivative_list)

        # ============================ online learning =============================

        initial_point = X_[np.random.randint(grid_size)]
        f_output.write("\n{0}, ".format(count))

        print ("\nGPUCB count:{0}...".format(count))
        GPUCBsolver = GPUCB(decomposition.get_function_value, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, delta=delta, a=a, b=b, X_=X_, initial_point=initial_point, discrete=discrete, linear=True, scale_down_factor=scale_down_factor, maxmin=maxmin) # linear arg only changes the beta_t used in exploration
        GPUCBsolver.run(total_run)
        GPUCB_scores[count] = GPUCBsolver.regret
        GPUCB_regret_list[count] = np.array(GPUCBsolver.regret_list)
        f_output.write("{0}, ".format(GPUCBsolver.regret))

        print ("\ndecomposed GPUCB count:{0}...".format(count))
        decomposedGPUCBsolver = DecomposedGPUCB(decomposition, kernelList, dimension, upper_bound, constraints, gp_alpha=gp_alpha_list, delta=delta, a=a, b=b, X_=X_, initial_point=initial_point, discrete=discrete, scale_down_factor=scale_down_factor, maxmin=maxmin)
        decomposedGPUCBsolver.run(total_run)
        decomposedGPUCB_scores[count] = decomposedGPUCBsolver.regret
        decomposedGPUCB_regret_list[count] = np.array(decomposedGPUCBsolver.regret_list)
        f_output.write("{0}, ".format(decomposedGPUCBsolver.regret))

        print ("\nExpected Improvement count:{0}...".format(count))
        EIsolver = Improvement(decomposition.get_function_value, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, method="EI", X_=X_, initial_point=initial_point, discrete=discrete, maxmin=maxmin)
        EIsolver.run(total_run)
        EI_scores[count] = EIsolver.regret
        EI_regret_list[count] = np.array(EIsolver.regret_list)
        f_output.write("{0}, ".format(EIsolver.regret))

        # print ("\ndecomposed EI count:{0}...".format(count))
        # decomposedEIsolver = DecomposedGPUCB(decomposition, kernelList, dimension, upper_bound, constraints, gp_alpha=gp_alpha_list, delta=delta, a=a, b=b, X_=X_, initial_point=initial_point, discrete=discrete, scale_down_factor=scale_down_factor, method='EI', maxmin=maxmin)
        # decomposedEIsolver.run(total_run)
        # decomposedEI_scores[count] = decomposedEIsolver.regret
        # decomposedEI_regret_list[count] = np.array(decomposedEIsolver.regret_list)
        # f_output.write("{0}, ".format(decomposedEIsolver.regret))

        print ("\nProbability of Improvement count:{0}...".format(count))
        POIsolver = Improvement(decomposition.get_function_value, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, method="POI", X_=X_, initial_point=initial_point, discrete=discrete, maxmin=maxmin)
        POIsolver.run(total_run)
        POI_scores[count] = POIsolver.regret
        POI_regret_list[count] = np.array(POIsolver.regret_list)
        f_output.write("{0}, ".format(POIsolver.regret))

        # print ("\ndecomposed POI count:{0}...".format(count))
        # decomposedPOIsolver = DecomposedGPUCB(decomposition, kernelList, dimension, upper_bound, constraints, gp_alpha=gp_alpha_list, delta=delta, a=a, b=b, X_=X_, initial_point=initial_point, discrete=discrete, scale_down_factor=scale_down_factor, method='POI', maxmin=maxmin)
        # decomposedPOIsolver.run(total_run)
        # decomposedPOI_scores[count] = decomposedPOIsolver.regret
        # decomposedPOI_regret_list[count] = np.array(decomposedPOIsolver.regret_list)
        # f_output.write("{0}, ".format(decomposedPOIsolver.regret))


    pickle.dump((GPUCB_regret_list, decomposedGPUCB_regret_list, EI_regret_list, decomposedEI_regret_list, POI_regret_list, decomposedPOI_regret_list), open(output_path+"regret_list_{0}.p".format(filename), 'wb'))

