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

    args = parser.parse_args()
    filename = args.name

    data_path = "flu/"
    output_path = "flu/result/"
    # ========================= experimental setting ========================
    J = 9
    budget = 0.2*J
    year_range = np.array([3, 3, 3, 10, 15, 15, 15, 5, 10])
    dimension = J
    lower_bound = 0.1
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
    max_derivative_list = 10

    # ================================ experimental design ===============================
    method = None
    optimize_kernel = True
    true_optimal = 0
    total_count = 1
    total_run = 200
    a_count = 1
    a_list = np.array([0.001]) * np.mean(max_derivative_list)
    #a_list = np.array([1e-5, 2e-5, 5e-5, 0.0001, 0.0002]) * np.mean(max_derivative_list)
    b_count = 1
    b_list = np.array([0.1])

    GPUCB_scores = np.zeros((a_count, b_count))
    decomposedGPUCB_scores = np.zeros((a_count, b_count))
    GPUCB_regret_list = np.zeros((a_count, b_count, total_run))
    decomposed_regret_list = np.zeros((a_count, b_count, total_run))


    def initial_point_generator():
        initial_point = np.random.rand(dimension)
        initial_point = initial_point / sum(weights * initial_point) * (budget - lower_bound * dimension) + lower_bound
        return initial_point


    for count in range(total_count):
        # GPUCB_scores = np.zeros((a_count, b_count))
        # decomposedGPUCB_scores = np.zeros((a_count, b_count))

        for a_index in range(a_count):
            a = a_list[a_index]
            for b_index in range(b_count):
                b = b_list[b_index]

                initial_point = initial_point_generator()

                print ("\nGPUCB count:{0}, a index:{1}, b index:{2}...".format(count, a_index, b_index))
                GPUCBsolver = GPUCB(decomposition.get_function_value, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, a=a, b=b, initial_point=initial_point, delta=delta, discrete=discrete, linear=linear, lower_bound=lower_bound, method=method, initial_point_generator=initial_point_generator, true_optimal=true_optimal, optimize_kernel=optimize_kernel) # linear arg only changes the beta_t used in exploration
                GPUCBsolver.run(total_run)
                GPUCB_scores[a_index, b_index] += GPUCBsolver.regret
                GPUCB_regret_list[a_index, b_index] += np.array(GPUCBsolver.regret_list)

                print ("\ndecomposed count:{0}, a index:{1}, b index:{2}...".format(count, a_index, b_index))
                decomposedGPUCBsolver = DecomposedGPUCB(decomposition, kernelList, dimension, upper_bound, constraints, gp_alpha=gp_alpha_list, a=a, b=b, initial_point=initial_point, delta=delta, discrete=discrete, lower_bound=lower_bound, method=method, initial_point_generator=initial_point_generator, true_optimal=true_optimal, optimize_kernel=optimize_kernel)
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

