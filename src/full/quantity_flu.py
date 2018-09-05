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
from flu import *

if __name__ == "__main__":
    # =========================== improvement quantity test ============================
    parser = argparse.ArgumentParser(description='Decomposed GPUCB and GPUCB comparison')
    parser.add_argument('-n', '--name', required=True, help='Input the name to save')
    parser.add_argument('-count', '--count', default=10, help='Input the total count')
    parser.add_argument('-iteration', '--iteration', required=True, help='Input the total iteration')

    args = parser.parse_args()
    filename = "{0}".format(args.name)
    total_count = int(args.count)
    total_iteration = int(args.iteration)

    data_path = "flu/"
    output_path = "quantity/flu/"

    # ========================= experimental setting ========================
    J = 5
    budget = 0.18 * 80
    #budget = 1.6
    #year_range = np.array([3, 3, 3, 10, 15, 15, 15, 5, 10])
    year_range = np.array([20, 30, 15, 5, 10])
    dimension = J
    lower_bound = 0.05
    upper_bound = 0.5
    gp_alpha_list = year_range * 0.0001 # TODO gp alpha list
    gp_alpha = sum(gp_alpha_list)
    delta = 0.1
    linear = True
    discrete = False

    # ============================== flu simulation setup ==============================
    transmissibility = 0.54
    contact_matrix = pd.read_csv(data_path+"contact.csv", index_col=0).values
    susceptibility = pd.read_csv(data_path+"susceptibility.csv").values.reshape(J)
    initial_population = pd.read_csv(data_path+"initial_population.csv").values.reshape(J)
    initial_infected = year_range * np.array([2, 1, 1.2, 1.5, 2]) * 91
    #initial_infected = [500, 300, 200, 800, 1200, 1000, 1500, 800, 2000]
    initial_recover = np.zeros(J)
    infectivity = np.ones(J) * 0.1

    beta_matrix = transmissibility * susceptibility.reshape((J,1)) * contact_matrix * infectivity
    death_rate = np.zeros(J)
    infected_death_rate = np.ones(J) * 0.00000008

    #weights = year_range * np.array([1, 1, 1, 1, 1, 1, 0.9, 0.6, 0.3])
    #weights = np.ones(J)
    weights = year_range * np.array([1.2, 1.1, 1, 0.8, 0.6])
    #weights = np.array([0.25, 0.40, 0.20, 0.1, 0.2])
    # constraints = [LinearConstraint([weights], [budget], [budget])]
    constraints = ({'type': 'eq', 'fun': lambda x: sum(x*weights) - budget})
    # constraints = None

    # plt.imshow(beta_matrix)
    # plt.show()

    # ================================= decomposition =====================================
    g = lambda x: sum(x)
    gList = [lambda x: 1 for i in range(J)]

    iterations = 3650
    decomposition = fluDecomposition(J, g, gList, beta_matrix, death_rate, infected_death_rate, initial_population, initial_infected, initial_recover, gp_alpha_list, iterations)

    sol = decomposition.simulation(np.ones(J)*0.2)

    S_list = np.sum(sol[:,:J], axis=1)
    I_list = np.sum(sol[:,J:J*2], axis=1)
    R_list = np.sum(sol[:,J*2:], axis=1)

    t_list = np.linspace(0, 365, iterations)

    def initial_point_generator():
        initial_point = np.random.rand(dimension)
        initial_point = initial_point / sum(weights * initial_point) * (budget - lower_bound * sum(weights)) + lower_bound
        return initial_point


    # ========================== experimental design ========================
    # ================================ experimental records ===============================

    f_output = open(output_path + "quantity_report_{0}.csv".format(filename), 'w')
    f_rmse = open(output_path + "quantity_rmse_{0}.csv".format(filename), 'w')
    for count in range(total_count):
        # ============================ generating kernels ==========================
        kernelList = []

        kernel_sample_size = 200
        print("roung {0}, fitting kernels...".format(count))
        X_ = np.zeros((kernel_sample_size, dimension))
        subfunction_values = np.zeros((kernel_sample_size, J))
        for i in range(kernel_sample_size):
            x = initial_point_generator()
            X_[i] = x
            subfunction_values[i] = decomposition.get_subfunction_values_no_noise(x)
        for i in range(J):
            # gpr = GaussianProcessRegressor(kernel=1.0*RBF(length_scale=1, length_scale_bounds=(2e-2, 2e2)), normalize_y=True)
            gpr = GaussianProcessRegressor(kernel=1.0 * Matern(length_scale=1, length_scale_bounds=(2e-2, 1), nu=1 + np.random.random()) +
                                                  #1.0 * Matern(length_scale=1, length_scale_bounds=(2e-2, 1), nu=1 + np.random.random()) + 
                                                  #1.0 * Matern(length_scale=1, length_scale_bounds=(2e-2, 1), nu=1 + np.random.random())
                                                  1.0 * RBF(length_scale=1, length_scale_bounds=(2e-2, 1)) +
                                                  1.0 * RationalQuadratic(alpha=0.5, length_scale=1, length_scale_bounds=(2e-2, 1))
                                                ,normalize_y=True)
            gpr.fit(X_, subfunction_values[:,i])
            kernelList.append(gpr.kernel_)

        kernel = sum(kernelList)
        print("kernel list:", kernelList)
        # ============================ generating functions ========================
        # random_seed = np.random.randint(10000)
        # print("random seed: {0}".format(random_seed))
        # targetList = np.array([GaussianProcessRegressor(kernel=kernelList[i], optimizer=None, alpha=gp_alpha).sample_y(X_, 1, random_state=random_seed) for i in range(J)])
        # fList = [randomizify(smoothify(targetList[i], X_), gp_alpha) for i in range(J)]

        # true_values = [g(targetList[:,i]) for i in range(grid_size)]

        # function_bounds = np.zeros(J)
        # for i in range(J):
        #     function_bounds[i] = np.mean(np.prod([np.abs(targetList[j]) for j in np.delete(np.arange(J), i)], axis=0))

        # b_scale = np.sum(function_bounds)
        # # gList = [lambda x: function_bounds[i] for i in range(J)]
        # gList = [lambda x: 1.0 for i in range(J)]

        # decomposition = Decomposition(J, fList, g, gList)

        # =============================== regression ===============================
        GPUCB_std_list = np.zeros(total_iteration)
        decomposedGPUCB_std_list = np.zeros(total_iteration)

        GPUCB_rmse_list = np.zeros(total_iteration)
        decomposedGPUCB_rmse_list = np.zeros(total_iteration)

        for iteration in range(total_iteration):
            print("count: {0}, iteration {1}".format(count, iteration))
            initial_point = [initial_point_generator() for i in range(iteration+1)]
            sample_sub_values = np.array([decomposition.get_subfunction_values(x) for x in initial_point])
            sample_value = np.array([decomposition.g(sub_value) for sub_value in sample_sub_values])

            whole_gpr = GaussianProcessRegressor(kernel=kernel, alpha=gp_alpha, optimizer=None, normalize_y=N_Y)
            whole_gpr.fit(initial_point, sample_value)

            gpr_list = []
            for i in range(J):
                tmp_gpr = GaussianProcessRegressor(kernel=kernelList[i], alpha=gp_alpha_list[i], optimizer=None, normalize_y=N_Y)
                tmp_gpr.fit(initial_point, sample_sub_values[:,i])
                gpr_list.append(tmp_gpr)

            random_size = 500
            random_points = [initial_point_generator() for i in range(random_size)] #[np.random.randint(0, grid_size, random_size)]
            true_values = [decomposition.get_function_value_no_noise(x) for x in random_points]

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
            print("GPUCB RMSE: {0}".format(GPUCB_rmse_list[iteration]))
            print("decomposed GPUCB RMSE: {0}".format(decomposedGPUCB_rmse_list[iteration]))

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


