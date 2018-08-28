import numpy as np
from scipy.integrate import odeint
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import argparse

from decomposition import *

# individual_normalization = np.array([0.46678064, 0.38247066, 0.30295984, 1.56013168, 1.67772638, 1.50348474, 2.29938285, 0.43154456, 1.08879793])
individual_normalization = np.array([0.4, 0.2, 0.2, 1, 1, 1, 1.5, 0.3, 0.7])


class fluDecomposition:
    def __init__(self, J, g, gList, beta_matrix, death_rate, infected_death_rate, initial_population, initial_infected, initial_recover, gp_alpha_list, iterations=3650):
        assert((J == beta_matrix.shape[0]) and (J == beta_matrix.shape[1]))
        self.J = J
        self.g = g
        self.gList = gList
        self.beta_matrix = beta_matrix
        self.death_rate = death_rate
        self.infected_death_rate = infected_death_rate
        self.initial_infected = initial_infected
        self.initial_population = initial_population
        self.initial_recover = initial_recover
        self.initial_susceptible = initial_population - initial_infected - initial_recover
        self.gp_alpha_list = gp_alpha_list
        self.iterations = iterations
        self.total_population = sum(initial_population)

    def derivative(self, y, t, v):
        assert(len(v) == self.J)
        beta = self.beta_matrix
        S = y[:self.J]
        I = y[self.J:self.J*2]
        R = y[self.J*2:]
        N = S + I + R
        I_proportion = I/N
        lamb = np.matmul(beta, I_proportion).reshape(self.J)
        dS = - lamb * S - self.death_rate * S
        dI = lamb * S - self.infected_death_rate * I - v * I
        dR = v * I - self.death_rate * R
        return np.concatenate([dS, dI, dR])

    def simulation(self, v):
        y0 = np.concatenate([self.initial_susceptible, self.initial_infected, self.initial_recover])
        t = np.linspace(0, 365, self.iterations) # iterations

        sol = odeint(self.derivative, y0, t, args=(v,))
        return sol

    def get_function_value(self, v):
        integration_I_list = self.get_subfunction_values(v)
        integration_I = g(integration_I_list)
        return integration_I

    def get_subfunction_values(self, v):
        sol = self.simulation(v)
        integration_I_list = np.sum(sol[:, self.J:self.J*2], axis=0) * 365 / (self.iterations * self.total_population) + (np.random.rand(self.J) - 0.5) * self.gp_alpha_list # - individual_normalization 
        return integration_I_list

    def get_coefficients(self, v):
        coefficients = np.zeros((self.J, 1))
        for i in range(self.J):
            coefficients[i] = self.gList[i](v)

        return coefficients


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Decomposed GPUCB and GPUCB comparison')
    parser.add_argument('-n', '--name', required=True, help='Input the name to save')
    parser.add_argument('-s', '--scale_down', default=5, help='Input the scale down factor')
    parser.add_argument('-iteration', '--iteration', default=300, help='Input the total iterations')
    parser.add_argument('-count', '--count', default=10, help='Input the total count')
    parser.add_argument('-a', '--a', required=True, help='Input the a value')
    parser.add_argument('-b', '--b', required=True, help='Input the b value')

    args = parser.parse_args()
    scale_down_factor = float(args.scale_down)
    filename = "{0}_scale{1}_a{2}_b{3}".format(args.name, args.scale_down, args.a, args.b)
    total_run = int(args.iteration)
    total_count = int(args.count)

    data_path = "flu/"
    output_path = "flu/new_result/"
    # ========================= experimental setting ========================
    J = 5
    budget = 0.15 * 80
    #budget = 1.6
    #year_range = np.array([3, 3, 3, 10, 15, 15, 15, 5, 10])
    year_range = np.array([20, 30, 15, 5, 10])
    dimension = J
    lower_bound = 0.05
    upper_bound = 0.6
    gp_alpha_list = year_range * 0.01 # TODO gp alpha list
    gp_alpha = sum(gp_alpha_list)
    delta = 0.05
    linear = True
    discrete = False

    # ============================== flu simulation setup ==============================
    transmissibility = 0.54
    contact_matrix = pd.read_csv(data_path+"contact.csv", index_col=0).values
    susceptibility = pd.read_csv(data_path+"susceptibility.csv").values.reshape(J)
    initial_population = pd.read_csv(data_path+"initial_population.csv").values.reshape(J)
    initial_infected = year_range * 91
    #initial_infected = [500, 300, 200, 800, 1200, 1000, 1500, 800, 2000]
    initial_recover = np.zeros(J)
    infectivity = np.ones(J) * 0.1

    beta_matrix = transmissibility * susceptibility.reshape((J,1)) * contact_matrix * infectivity
    death_rate = np.zeros(J)
    infected_death_rate = np.ones(J) * 0.0000008

    #weights = year_range * np.array([0.5, 0.7, 1, 1.2, 1.5, 1.8, 1.9, 2, 2.5])
    #weights = np.ones(J)
    weights = year_range * np.array([0.6, 0.7, 1, 1.2, 1.3])
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

    # plt.plot(t_list, S_list, 'b')
    # plt.plot(t_list, I_list, 'r')
    # plt.plot(t_list, R_list, 'g')
    # plt.show()


    # =================================== kernels determination =========================

    def initial_point_generator():
        initial_point = np.random.rand(dimension)
        initial_point = initial_point / sum(weights * initial_point) * (budget - lower_bound * sum(weights)) + lower_bound
        return initial_point

    #"""
    kernelList = []

    kernel_sample_size = 1000
    X_ = np.zeros((kernel_sample_size, dimension))
    subfunction_values = np.zeros((kernel_sample_size, J))
    for i in range(kernel_sample_size):
        x = initial_point_generator()
        X_[i] = x
        subfunction_values[i] = decomposition.get_subfunction_values(x)
    for i in range(J):
        # gpr = GaussianProcessRegressor(kernel=1.0*Matern(length_scale=1, length_scale_bounds=(2e-2, 2e2)), normalize_y=True)
        gpr = GaussianProcessRegressor(kernel=1.0*RBF(length_scale=1, length_scale_bounds=(2e-2, 2e2)) + 1.0*Matern(length_scale=1, length_scale_bounds=(2e-2, 2e2)) + 1.0 *
                + 1.0*RationalQuadratic(alpha=0.1, length_scale=1, length_scale_bounds=(2e-2, 1)), normalize_y=True)
        gpr.fit(X_, subfunction_values[:,i])
        kernelList.append(gpr.kernel_)

    print("kernel list: {0}".format(kernelList))
    #"""

    # ==================================== kernels =======================================

    # kernelList = [1/float(J) * RBF(length_scale=1, length_scale_bounds=(2e-2, 2e1)) for i in range(J)]
    # kernel_variance = np.array([ 8.23**2, 10.7**2, 3.69**2, 0.699**2, 2.08**2])
    # kernelList = np.array([RBF(length_scale=0.399), RBF(length_scale=0.434), RBF(length_scale=0.347),
    #                        RBF(length_scale=0.431), RBF(length_scale=0.4), ]) * kernel_variance

    # kernel_variance = np.array([ 2.53**2, 1.92**2, 2.27**2, 4**2, 5.98**2, 5.47**2, 3.92**2, 0.903**2, 3.43**2 ])
    # kernelList = np.array([RBF(length_scale=1.790), RBF(length_scale=1.760), RBF(length_scale=1.860),
    #                       RBF(length_scale=1), RBF(length_scale=0.697), RBF(length_scale=0.680),
    #                       RBF(length_scale=0.770), RBF(length_scale=1.510), RBF(length_scale=2.280), ]) * kernel_variance
    # kernel = RBF(length_scale=0.5, length_scale_bounds=(1e-1, 1e1))
    kernel = sum(kernelList)

    max_derivative_list = 10

    # ================================ experimental design ===============================
    optimization_method = None
    optimize_kernel = False
    true_optimal = 0 # Empirically observed
    a = float(args.a)
    b = float(args.b) * max_derivative_list

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
        # GPUCB_scores = np.zeros((a_count, b_count))
        # decomposedGPUCB_scores = np.zeros((a_count, b_count))

        initial_point = initial_point_generator()
        f_output.write("\n{0}, ".format(count))

        print ("\nGPUCB count: {0}...".format(count))
        GPUCBsolver = GPUCB(decomposition.get_function_value, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, a=a, b=b, initial_point=initial_point, delta=delta, discrete=discrete, linear=linear, lower_bound=lower_bound, optimization_method=optimization_method, initial_point_generator=initial_point_generator, true_optimal=true_optimal, optimize_kernel=optimize_kernel, scale_down_factor=scale_down_factor, maxmin="min") # linear arg only changes the beta_t used in exploration
        GPUCBsolver.run(total_run)
        GPUCB_scores[count] = GPUCBsolver.regret
        GPUCB_regret_list[count] = np.array(GPUCBsolver.regret_list)
        f_output.write("{0}, ".format(GPUCBsolver.regret))

        print ("\ndecomposed GPUCB count: {0}...".format(count))
        decomposedGPUCBsolver = DecomposedGPUCB(decomposition, kernelList, dimension, upper_bound, constraints, gp_alpha=gp_alpha_list, a=a, b=b, initial_point=initial_point, delta=delta, discrete=discrete, lower_bound=lower_bound, optimization_method=optimization_method, initial_point_generator=initial_point_generator, true_optimal=true_optimal, optimize_kernel=optimize_kernel, scale_down_factor=scale_down_factor, maxmin="min")
        decomposedGPUCBsolver.run(total_run)
        decomposedGPUCB_scores[count] = decomposedGPUCBsolver.regret
        decomposedGPUCB_regret_list[count] = np.array(decomposedGPUCBsolver.regret_list)
        f_output.write("{0}, ".format(decomposedGPUCBsolver.regret))

        print ("\nExpected Improvement count:{0}...".format(count))
        EIsolver = Improvement(decomposition.get_function_value, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, method="EI", initial_point=initial_point, discrete=discrete, lower_bound=lower_bound, optimization_method=optimization_method, initial_point_generator=initial_point_generator, true_optimal=true_optimal, optimize_kernel=optimize_kernel, maxmin="min")
        EIsolver.run(total_run)
        EI_scores[count] = EIsolver.regret
        EI_regret_list[count] = np.array(EIsolver.regret_list)
        f_output.write("{0}, ".format(EIsolver.regret))

        # print ("\ndecomposed EI count: {0}...".format(count))
        # decomposedEIsolver = DecomposedGPUCB(decomposition, kernelList, dimension, upper_bound, constraints, gp_alpha=gp_alpha_list, a=a, b=b, initial_point=initial_point, delta=delta, discrete=discrete, lower_bound=lower_bound, optimization_method=optimization_method, initial_point_generator=initial_point_generator, true_optimal=true_optimal, optimize_kernel=optimize_kernel, method="EI", maxmin="min")
        # decomposedEIsolver.run(total_run)
        # decomposedEI_scores[count] = decomposedEIsolver.regret
        # decomposedEI_regret_list[count] = np.array(decomposedEIsolver.regret_list)
        # f_output.write("{0}, ".format(decomposedEIsolver.regret))

        print ("\nProbability of Improvement count:{0}...".format(count))
        POIsolver = Improvement(decomposition.get_function_value, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, method="POI", initial_point=initial_point, discrete=discrete, lower_bound=lower_bound, optimization_method=optimization_method, initial_point_generator=initial_point_generator, true_optimal=true_optimal, optimize_kernel=optimize_kernel, maxmin="min")
        POIsolver.run(total_run)
        POI_scores[count] = POIsolver.regret
        POI_regret_list[count] = np.array(POIsolver.regret_list)
        f_output.write("{0}, ".format(POIsolver.regret))

        # print ("\ndecomposed POI count: {0}...".format(count))
        # decomposedPOIsolver = DecomposedGPUCB(decomposition, kernelList, dimension, upper_bound, constraints, gp_alpha=gp_alpha_list, a=a, b=b, initial_point=initial_point, delta=delta, discrete=discrete, lower_bound=lower_bound, optimization_method=optimization_method, initial_point_generator=initial_point_generator, true_optimal=true_optimal, optimize_kernel=optimize_kernel, method="POI", maxmin="min")
        # decomposedPOIsolver.run(total_run)
        # decomposedPOI_scores[count] = decomposedPOIsolver.regret
        # decomposedPOI_regret_list[count] = np.array(decomposedPOIsolver.regret_list)
        # f_output.write("{0}".format(decomposedPOIsolver.regret))

        pickle.dump((GPUCB_regret_list, decomposedGPUCB_regret_list, EI_regret_list, decomposedEI_regret_list, POI_regret_list, decomposedPOI_regret_list), open(output_path+"regret_list_{0}.p".format(filename), 'wb'))

    f_output.close()


