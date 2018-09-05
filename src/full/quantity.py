import numpy as np
import scipy
from sklearn.gaussian_process.kernels import (RBF, Matern, RationalQuadratic,
                                              ExpSineSquared, DotProduct,
                                              ConstantKernel, WhiteKernel,
                                              RationalQuadratic)

from sklearn.metrics import mean_squared_error

import pandas as pd
import argparse
import pickle

from decomposition import *

if __name__ == "__main__":
    # =========================== improvement quantity test ============================
    parser = argparse.ArgumentParser(description='Decomposed GPUCB and GPUCB comparison')
    parser.add_argument('-n', '--name', required=True, help='Input the name to save')
    parser.add_argument('-J', '--J', required=True, help='Input the number of decompositions')
    parser.add_argument('-count', '--count', default=10, help='Input the total count')
    parser.add_argument('-iteration', '--iteration', required=True, help='Input the total iteration')
    parser.add_argument('-kernel', '--kernel', required=True, help='Input the kernel')

    args = parser.parse_args()
    filename = "{0}_J{1}".format(args.name, args.J)
    total_count = int(args.count)
    total_iteration = int(args.iteration)

    output_path = "quantity/result/"

    J = int(args.J)
    upper_bound = 1
    grid_size = 2000
    gp_alpha_list = [0.0001] * J
    gp_alpha = 0.0001 * J
    kernel_name = args.kernel

    X_ = np.reshape(np.linspace(0, upper_bound, grid_size), (grid_size, 1))

    g = lambda x: np.sum(x)

    # ========================== experimental design ========================
    # ================================ experimental records ===============================

    f_output = open(output_path + "quantity_report_{0}.csv".format(filename), 'w')
    f_rmse = open(output_path + "quantity_rmse_{0}.csv".format(filename), 'w')
    for count in range(total_count):
        # ============================ generating kernels ==========================
        if kernel_name == "RQ":
            kernelList = [RationalQuadratic(alpha=0.05 + 0.2 * np.random.random(), length_scale=0.5 + 1.5*np.random.random()) for i in range(1,J+1)]
        elif kernel_name == "Matern":
            kernelList = [Matern(length_scale=0.5 + 1.5*np.random.random(), nu=1 + np.random.random()) for i in range(1,J+1)]
        elif kernel_name == "RBF":
            kernelList = [RBF(length_scale=0.5 + 1.5*np.random.random()) for i in range(1,J+1)]
        else:
            raise(Exception("Not supported kernel!"))

        kernel = sum(kernelList)
        print("kernel list:", kernelList)
        # ============================ generating functions ========================
        random_seed = np.random.randint(10000)
        print("random seed: {0}".format(random_seed))
        targetList = np.array([GaussianProcessRegressor(kernel=kernelList[i], optimizer=None, alpha=gp_alpha).sample_y(X_, 1, random_state=random_seed) for i in range(J)])
        fList = [randomizify(smoothify(targetList[i], X_), gp_alpha) for i in range(J)]

        true_values = [g(targetList[:,i]) for i in range(grid_size)]

        function_bounds = np.zeros(J)
        for i in range(J):
            function_bounds[i] = np.mean(np.prod([np.abs(targetList[j]) for j in np.delete(np.arange(J), i)], axis=0))

        b_scale = np.sum(function_bounds)
        # gList = [lambda x: function_bounds[i] for i in range(J)]
        gList = [lambda x: 1.0 for i in range(J)]

        decomposition = Decomposition(J, fList, g, gList)

        # =============================== regression ===============================
        GPUCB_std_list = np.zeros(total_iteration)
        decomposedGPUCB_std_list = np.zeros(total_iteration)

        GPUCB_rmse_list = np.zeros(total_iteration)
        decomposedGPUCB_rmse_list = np.zeros(total_iteration)

        for iteration in range(total_iteration):
            print("count: {0}, iteration {1}".format(count, iteration))
            initial_point = X_[np.random.randint(0, grid_size, iteration+1)]
            sample_sub_values = np.array([decomposition.get_subfunction_values(x) for x in initial_point])
            sample_value = np.array([decomposition.g(sub_value) for sub_value in sample_sub_values])

            whole_gpr = GaussianProcessRegressor(kernel=kernel, alpha=gp_alpha, optimizer=None, normalize_y=N_Y)
            whole_gpr.fit(initial_point, sample_value)

            gpr_list = []
            for i in range(J):
                tmp_gpr = GaussianProcessRegressor(kernel=kernelList[i], alpha=gp_alpha_list[i], optimizer=None, normalize_y=N_Y)
                tmp_gpr.fit(initial_point, sample_sub_values[:,i])
                gpr_list.append(tmp_gpr)

            random_size = grid_size
            random_points = X_ #[np.random.randint(0, grid_size, random_size)]

            whole_mean, whole_std = whole_gpr.predict(random_points, return_std=True)
            decomposed_mean = np.zeros(random_size)
            decomposed_variance = np.zeros(random_size)
            for i in range(J):
                tmp_mean, tmp_std = gpr_list[i].predict(random_points, return_std=True)
                decomposed_mean += tmp_mean
                decomposed_variance += np.square(tmp_std)
            decomposed_std = np.sqrt(decomposed_variance)

            GPUCB_std_list[iteration] = np.mean(whole_std)
            decomposedGPUCB_std_list[iteration] = np.mean(decomposed_std)

            GPUCB_rmse_list[iteration] = mean_squared_error(whole_mean, true_values)
            decomposedGPUCB_rmse_list[iteration] = mean_squared_error(decomposed_mean, true_values)

        improvement_ratio = GPUCB_std_list / decomposedGPUCB_std_list
        improvement_rmse_ratio = GPUCB_rmse_list / decomposedGPUCB_rmse_list

        f_output.write("GPUCB, count, {0},".format(count) + ",".join([str(x) for x in GPUCB_std_list]) + "\n")
        f_output.write("decomposed GPUCB, count, {0},".format(count) + ",".join([str(x) for x in decomposedGPUCB_std_list]) + "\n")
        f_output.write("improvement ratio, count, {0},".format(count) + ",".join([str(x) for x in improvement_ratio]) + "\n")

        f_rmse.write("GPUCB, count, {0},".format(count) + ",".join([str(x) for x in GPUCB_rmse_list]) + "\n")
        f_rmse.write("decomposed GPUCB, count, {0},".format(count) + ",".join([str(x) for x in decomposedGPUCB_rmse_list]) + "\n")
        f_rmse.write("improvement ratio, count, {0},".format(count) + ",".join([str(x) for x in improvement_rmse_ratio]) + "\n")


    f_output.close()
    f_rmse.close()


