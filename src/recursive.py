import matlab.engine
import os
import uptakePolicy
import numpy as np
import itertools
import pickle
import sys
import argparse
import cplex
import derivative
import util

import params

# =============== matplot ================
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

# ============= our library ==============
Bryan_path = "/home/kai/Dropbox/USC/publication/TBsimulationCode/Bryan/"
dir_path = "/home/kai/Dropbox/USC/publication/TBsimulationCode/"
sys.path.insert(0, Bryan_path)

import sfw

# ============== global variables =============
MAX_AGE = 110

if __name__ == "__main__":

    yrsOfAnalysis = 25
    n = 110
    effort = 10
    random_iterations = 20
    recursive_iterations = 1
    restart = False

    newE, newI, alpha_fast, alpha_slow, G, mu, death_rate = util.get_parameters(n, yrsOfAnalysis)

    parser = argparse.ArgumentParser(description="filename parser")
    parser.add_argument("-d", "--date", help="Input date")
    parser.add_argument("-m", "--method", help="Input one of the methods: initial, optimal, full, peak40, peak55")
    parser.add_argument("-n", "--number", help="Input the number of initial population")
    parser.add_argument("-v", "--variance", default=0, help="Input the change to the rerun nu")
    parser.add_argument("-t", "--testingSample", default=0, help="Input the number of testing sample")

    args = parser.parse_args()
    print "args: {0}".format(args)

    date = args.date
    file_index = args.method
    NumPpl = int(args.number)
    rerun_variance = float(args.variance)
    testing_sample = int(args.testingSample)

    print "date: {0}, number of population {1}".format(date, NumPpl)

    # ============================= recording ==============================
    low_objective_value_list = []
    high_objective_value_list = []

    # ====================== loading previous data =========================
    print "======================== Loading ... ============================"
    if not restart:
        normalized_high_infected_list, normalized_high_population_list, nu_list = pickle.load(open(dir_path + "data/high_fidelity_{0}_{1}_{2}.data".format(date, file_index, NumPpl), "rb"))
    else:
        normalized_high_infected_list = []
        normalized_high_population_list = []
        nu_list = []
        for i in range(random_iterations):
            nu = util.get_random_nu(effort, n)

            tmp_nu = np.resize(nu, (n,1))
            tmp_high_infected, tmp_high_population = util.query_high(NumPpl, 25, np.concatenate([tmp_nu, tmp_nu], axis=1).tolist(), uptakePolicy.finerUptakeAgeBracs)

            tmp_high_infected = np.array(tmp_high_infected)
            tmp_high_population = np.array(tmp_high_population)

            normalized_high_infected = tmp_high_infected[:, range(130*12, 155*12, 12)]
            normalized_high_population = tmp_high_population[:, range(130*12, 155*12, 12)]

            normalized_high_infected_list.append(normalized_high_infected)
            normalized_high_population_list.append(normalized_high_population)
            nu_list.append(nu)


    # =========================== recursive running ============================
    for iteration in range(recursive_iterations):
        # ======================== initial population ==========================
        I_infected = np.zeros((n,n))
        N_population = np.zeros((n,n))
        S_safe = np.zeros((n,n))

        initial_population = np.zeros(n)
        initial_infected = np.zeros(n)
        previous_sample_size = len(normalized_high_infected_list)
        for i in range(previous_sample_size):
            initial_population += normalized_high_population_list[i][:,0]
            initial_infected += normalized_high_infected_list[i][:,0]
        initial_population /= previous_sample_size
        initial_infected /= previous_sample_size

        I_infected[:,0] = initial_infected
        N_population[:,0] = initial_population
        S_safe[:,0] = initial_population - initial_infected

        # ===================== extract the constraints ========================
        print "=================== extracting constraints ======================"
        x_train_list = []
        y_train_list = []

        samples_number = len(normalized_high_infected_list)
        #samples_number = 1
        birth_rate = sfw.zeros((n, yrsOfAnalysis))
        newI = np.zeros((n, yrsOfAnalysis))
        for k in range(samples_number):
            tmp_high_infected = normalized_high_infected_list[k]
            tmp_high_population = normalized_high_population_list[k]
            tmp_birth_rate = sfw.zeros((n, yrsOfAnalysis))
            tmp_newI = np.zeros((n, yrsOfAnalysis))
            for i in range(yrsOfAnalysis-1):
                tmp_birth_rate[0, i] = tmp_high_population[0,i+1]
                tmp_newI[0,i] = tmp_high_infected[0,i+1]
            birth_rate += tmp_birth_rate
            newI += tmp_newI
            #new_x, new_y = find_beta_constraints(yrsOfAnalysis, tmp_high_infected, tmp_high_population, tmp_newI, mu, nu_list[k], death_rate, tmp_birth_rate) # remember to modify "normalized..."
            new_x, new_y = util.separated_find_beta_constraints(yrsOfAnalysis, tmp_high_infected, tmp_high_population, tmp_newI, mu, nu_list[k], death_rate, tmp_birth_rate) # remember to modify "normalized..."

            x_train_list.append(new_x)
            y_train_list.append(new_y)

        birth_rate /= float(samples_number)
        newI /= float(samples_number)

        # ============================= training =============================
        print "======================== training ============================="
        x_train = np.concatenate(x_train_list, axis=1)
        y_train = np.concatenate(y_train_list, axis=1)

        weight = np.ones((n, samples_number * (yrsOfAnalysis - 1))) # TODO
        #for i in range(n):
        #    for j in range(samples_number * (yrsOfAnalysis - 1)):
        #        weight[i][j] = float(high_population[i][j % (yrsOfAnalysis - 1)] ** 2)

        print "====================== running linear regression ======================="
        beta = util.quadratic_solver(x_train, y_train, weight)

        # =============== optimization of low fidelity model =================

        L = np.ones(110) * 0.01
        U = np.ones(110) * 0.5
        K = effort - np.sum(L)
        num_iter = 10
        c = np.bmat([[np.ones((n, 1))], [np.zeros((n+1, 1))]])

        nu = np.zeros((n))
        for i in range(num_iter):
            print i

            x0 = np.transpose(np.bmat([np.transpose(I_infected[:,0]), np.zeros((n)), np.ones((1))]))

#            grad = gradient(n, T, G, S[j], x0, c, newE, newI, nu+L, mu, d, alpha_fast, alpha_slow, beta[j], N[j])
            #Swap this line in for Tth step only, vs sum of 1...T 
#            grad = grad + gradient(n, T, G, S[j], x0, c, newE, newI, nu+L, mu, d, alpha_fast, alpha_slow, beta[j], N[j])

            grad = sfw.gradient_sum(n, yrsOfAnalysis, G, S_safe, x0, c, newE, newI, nu+L, mu, death_rate, alpha_fast, alpha_slow, beta, N_population)

            #derI, derS = derivative.getDerivativeByNu(n, yrsOfAnalysis, I_infected, S_safe, newI, birth_rate, death_rate, nu, mu, beta)
            #grad = - np.sum(derI, axis=(1,2))

            #grad = grad / len(beta)
            v = sfw.greedy(grad, U - L, np.zeros((n)), K)
            nu = nu + (1./num_iter)*v

            low_infected, low_population = util.query_low(n, yrsOfAnalysis, G, S_safe, newE, newI, nu, mu, death_rate, alpha_fast, alpha_slow, beta, N_population, I_infected, birth_rate)
            tmp_objective_value = np.sum(low_infected)
            print "objective value at round {0}: {1}".format(i, tmp_objective_value)

        nu = nu + L

        # ================ running both fidelity models to check =============
        # ----------------------- high fidelity model ------------------------
        """
        tmp_nu = np.resize(nu, (n,1))
        tmp_high_infected, tmp_high_population = util.query_high(NumPpl, 25, np.concatenate([tmp_nu, tmp_nu], axis=1).tolist(), uptakePolicy.finerUptakeAgeBracs)

        tmp_high_infected = np.array(tmp_high_infected)
        tmp_high_population = np.array(tmp_high_population)

        normalized_high_infected = tmp_high_infected[:, range(130*12, 155*12, 12)]
        normalized_high_population = tmp_high_population[:, range(130*12, 155*12, 12)]

        normalized_high_infected_list.append(normalized_high_infected)
        normalized_high_population_list.append(normalized_high_population)

        nu_list.append(nu)

        high_objective_value = np.sum(normalized_high_infected) 

        # ----------------------- low fidelity model -------------------------
        low_infected, low_population = util.query_low(n, yrsOfAnalysis, G, S_safe, newE, newI, nu, mu, death_rate, alpha_fast, alpha_slow, beta, N_population, I_infected, birth_rate) # query_low(n, optTimeHorizon, G, S, newE, newI, nu, mu, d, alpha_fast, alpha_slow, beta, N, I, b):
        low_objective_value = np.sum(low_infected)

        # ----------------------------- recording ----------------------------
        low_objective_value_list.append(low_objective_value)
        high_objective_value_list.append(high_objective_value)

        print "iteration: {0}, high objective value: {1}, low objective value: {2}".format(iteration, high_objective_value, low_objective_value)
        infected_variation = np.sum(np.abs(normalized_high_infected - low_infected))
        population_variation = np.sum(np.abs(normalized_high_population - low_population))
        print "total variation infected: {0}, population: {1}".format(infected_variation, population_variation)
        

        # ============================= pickle save ==============================
        pickle.dump((normalized_high_infected_list, normalized_high_population_list, nu_list), open(dir_path + "data/high_fidelity_{0}_{1}_{2}.data".format(date, file_index, NumPpl), "wb"))
        """
