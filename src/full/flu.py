import numpy as np
from scipy.integrate import odeint
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import argparse
from scipy.optimize import LinearConstraint

from decomposition import *


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
        sol = self.simulation(v)
        integration_I_list = np.sum(sol[:, self.J:self.J*2], axis=0) * 365 / (self.iterations * self.total_population) + np.random.rand(self.J) * self.gp_alpha_list
        integration_I = g(integration_I_list)
        return -integration_I

    def get_subfunction_values(self, v):
        sol = self.simulation(v)
        integration_I_list = np.sum(sol[:, self.J:self.J*2], axis=0) * 365 / (self.iterations * self.total_population) + np.random.rand(self.J) * self.gp_alpha_list
        return -integration_I_list

    def get_coefficients(self, v):
        coefficients = np.zeros((self.J, 1))
        for i in range(self.J):
            coefficients[i] = self.gList[i](v)

        return coefficients


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Decomposed GPUCB and GPUCB comparison')
    parser.add_argument('-n', '--name', required=True, help='Input the name to save')
    parser.add_argument('-s', '--scale_down', required=True, help='Input the scale down factor')
    parser.add_argument('-iteration', '--iteration', default=300, help='Input the total iterations')
    parser.add_argument('-count', '--count', default=10, help='Input the total count')

    args = parser.parse_args()
    scale_down_factor = float(args.scale_down)
    filename = "{0}_scale{1}".format(args.name, args.scale_down)
    total_run = int(args.iteration)
    total_count = int(args.count)

    data_path = "flu/"
    output_path = "flu/new_result/"
    # ========================= experimental setting ========================
    J = 9
    budget = 0.2*J
    year_range = np.array([3, 3, 3, 10, 15, 15, 15, 5, 10])
    dimension = J
    lower_bound = 0.05
    upper_bound = 1
    weights = np.ones(J)
    # constraints = [LinearConstraint([weights], [budget], [budget])]
    constraints = ({'type': 'eq', 'fun': lambda x: sum(x*weights) - budget})
    # constraints = None
    gp_alpha_list = year_range * 0.001 # TODO gp alpha list
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
    initial_recover = np.zeros(J)
    infectivity = np.ones(J) * 0.1

    beta_matrix = transmissibility * susceptibility.reshape((J,1)) * contact_matrix * infectivity
    death_rate = np.zeros(J)
    infected_death_rate = np.ones(J) * 0.000016

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

    # ==================================== kernels =======================================

    kernelList = [RBF(length_scale=0.5, length_scale_bounds=(1e-1, 1e1)) for i in range(J)]
    kernel = RBF(length_scale=0.5, length_scale_bounds=(1e-1, 1e1))

    function_bounds = np.ones(J)
    max_derivative_list = 1

    # ================================ experimental design ===============================
    optimization_method = None
    optimize_kernel = True
    true_optimal = -2.45 # Empirically observed
    a = np.mean(max_derivative_list)
    b = 1

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

    def initial_point_generator():
        initial_point = np.random.rand(dimension)
        initial_point = initial_point / sum(weights * initial_point) * (budget - lower_bound * dimension) + lower_bound
        return initial_point


    f_output = open(output_path + "scores_report_{0}.csv".format(filename), 'w')
    f_output.write("round, GPUCB, decomposed GPUCB, EI, decomposed EI, POI, decomposed POI")
    for count in range(total_count):
        # GPUCB_scores = np.zeros((a_count, b_count))
        # decomposedGPUCB_scores = np.zeros((a_count, b_count))

        initial_point = initial_point_generator()
        f_output.write("\n{0}, ".format(count))

        print ("\nGPUCB count: {0}...".format(count))
        GPUCBsolver = GPUCB(decomposition.get_function_value, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, a=a, b=b, initial_point=initial_point, delta=delta, discrete=discrete, linear=linear, lower_bound=lower_bound, optimization_method=optimization_method, initial_point_generator=initial_point_generator, true_optimal=true_optimal, optimize_kernel=optimize_kernel, scale_down_factor=scale_down_factor) # linear arg only changes the beta_t used in exploration
        GPUCBsolver.run(total_run)
        GPUCB_scores[count] = GPUCBsolver.regret
        GPUCB_regret_list[count] = np.array(GPUCBsolver.regret_list)
        f_output.write("{0}, ".format(GPUCBsolver.regret))

        print ("\ndecomposed GPUCB count: {0}...".format(count))
        decomposedGPUCBsolver = DecomposedGPUCB(decomposition, kernelList, dimension, upper_bound, constraints, gp_alpha=gp_alpha_list, a=a, b=b, initial_point=initial_point, delta=delta, discrete=discrete, lower_bound=lower_bound, optimization_method=optimization_method, initial_point_generator=initial_point_generator, true_optimal=true_optimal, optimize_kernel=optimize_kernel, scale_down_factor=scale_down_factor)
        decomposedGPUCBsolver.run(total_run)
        decomposedGPUCB_scores[count] = decomposedGPUCBsolver.regret
        decomposedGPUCB_regret_list[count] = np.array(decomposedGPUCBsolver.regret_list)
        f_output.write("{0}, ".format(decomposedGPUCBsolver.regret))

        print ("\nExpected Improvement count:{0}...".format(count))
        EIsolver = Improvement(decomposition.get_function_value, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, method="EI", initial_point=initial_point, discrete=discrete, lower_bound=lower_bound, optimization_method=optimization_method, initial_point_generator=initial_point_generator, true_optimal=true_optimal, optimize_kernel=optimize_kernel)
        EIsolver.run(total_run)
        EI_scores[count] = EIsolver.regret
        EI_regret_list[count] = np.array(EIsolver.regret_list)
        f_output.write("{0}, ".format(EIsolver.regret))

        print ("\ndecomposed EI count: {0}...".format(count))
        decomposedEIsolver = DecomposedGPUCB(decomposition, kernelList, dimension, upper_bound, constraints, gp_alpha=gp_alpha_list, a=a, b=b, initial_point=initial_point, delta=delta, discrete=discrete, lower_bound=lower_bound, optimization_method=optimization_method, initial_point_generator=initial_point_generator, true_optimal=true_optimal, optimize_kernel=optimize_kernel, method="EI")
        decomposedEIsolver.run(total_run)
        decomposedEI_scores[count] = decomposedEIsolver.regret
        decomposedEI_regret_list[count] = np.array(decomposedEIsolver.regret_list)
        f_output.write("{0}, ".format(decomposedEIsolver.regret))

        print ("\nProbability of Improvement count:{0}...".format(count))
        POIsolver = Improvement(decomposition.get_function_value, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, method="POI", initial_point=initial_point, discrete=discrete, lower_bound=lower_bound, optimization_method=optimization_method, initial_point_generator=initial_point_generator, true_optimal=true_optimal, optimize_kernel=optimize_kernel)
        POIsolver.run(total_run)
        POI_scores[count] = POIsolver.regret
        POI_regret_list[count] = np.array(POIsolver.regret_list)
        f_output.write("{0}, ".format(POIsolver.regret))

        print ("\ndecomposed POI count: {0}...".format(count))
        decomposedPOIsolver = DecomposedGPUCB(decomposition, kernelList, dimension, upper_bound, constraints, gp_alpha=gp_alpha_list, a=a, b=b, initial_point=initial_point, delta=delta, discrete=discrete, lower_bound=lower_bound, optimization_method=optimization_method, initial_point_generator=initial_point_generator, true_optimal=true_optimal, optimize_kernel=optimize_kernel, method="POI")
        decomposedPOIsolver.run(total_run)
        decomposedPOI_scores[count] = decomposedPOIsolver.regret
        decomposedPOI_regret_list[count] = np.array(decomposedPOIsolver.regret_list)
        f_output.write("{0}".format(decomposedPOIsolver.regret))

    f_output.close()

    pickle.dump((GPUCB_regret_list, decomposedGPUCB_regret_list, EI_regret_list, decomposedEI_regret_list, POI_regret_list, decomposedPOI_regret_list), open(output_path+"regret_list_{0}.p".format(filename), 'wb'))

