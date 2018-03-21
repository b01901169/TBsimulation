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
    random_iterations = 5
    recursive_iterations = 100
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
        normalized_high_infected_list, normalized_high_latent_list, normalized_high_population_list, nu_list = pickle.load(open(dir_path + "data/high_fidelity_{0}_{1}_{2}.data".format(date, file_index, NumPpl), "rb"))
    else:
        normalized_high_infected_list = []
        normalized_high_latent_list = []
        normalized_high_population_list = []
        nu_list = []
        for i in range(random_iterations):
            nu = util.get_random_nu(effort, n)
            if i == random_iterations - 1:
                nu = np.ones(n) * effort / n

            tmp_nu = np.resize(nu, (n,1))
            tmp_high_infected, tmp_high_latent, tmp_high_population = util.query_high(NumPpl, 25, np.concatenate([tmp_nu, tmp_nu], axis=1).tolist(), uptakePolicy.finerUptakeAgeBracs)

            tmp_high_infected = np.array(tmp_high_infected)
            tmp_high_latent = np.array(tmp_high_latent)
            tmp_high_population = np.array(tmp_high_population)

            normalized_high_infected = tmp_high_infected[:, range(130*12, 155*12, 12)]
            normalized_high_latent = tmp_high_latent[:, range(130*12, 155*12, 12)]
            normalized_high_population = tmp_high_population[:, range(130*12, 155*12, 12)]

            normalized_high_infected_list.append(normalized_high_infected)
            normalized_high_latent_list.append(normalized_high_latent)
            normalized_high_population_list.append(normalized_high_population)
            nu_list.append(nu)

        pickle.dump((normalized_high_infected_list, normalized_high_latent_list, normalized_high_population_list, nu_list), open(dir_path + "data/high_fidelity_{0}_{1}_{2}.data".format(date, file_index, NumPpl), "wb"))


    # =========================== recursive running ============================
    for iteration in range(recursive_iterations):
        # ======================== initial population ==========================
        I_infected = np.zeros((n,yrsOfAnalysis))
        E_latent = np.zeros((n,yrsOfAnalysis))
        N_population = np.zeros((n,yrsOfAnalysis))
        S_safe = np.zeros((n,yrsOfAnalysis))

        initial_population = np.zeros(n)
        initial_latent = np.zeros(n)
        initial_infected = np.zeros(n)
        previous_sample_size = len(normalized_high_infected_list)
        for i in range(previous_sample_size):
            initial_population += normalized_high_population_list[i][:,0]
            initial_latent += normalized_high_latent_list[i][:,0]
            initial_infected += normalized_high_infected_list[i][:,0]
        initial_population /= previous_sample_size
        initial_latent /= previous_sample_size
        initial_infected /= previous_sample_size

        I_infected[:,0] = initial_infected
        E_latent[:,0] = initial_latent
        N_population[:,0] = initial_population
        S_safe[:,0] = initial_population - initial_infected - initial_latent

        # ===================== extract the constraints ========================
        print "=================== extracting constraints ======================"
        x_train_infected_list = []
        y_train_infected_list = []
        x_train_latent_list = []
        y_train_latent_list = []

        samples_number = len(normalized_high_infected_list)
        train_set = range(max(0,samples_number-10),samples_number) # TODO
        train_number = len(train_set)
        #samples_number = 1
        birth_rate = sfw.zeros((n, yrsOfAnalysis))
        newI = np.zeros((n, yrsOfAnalysis))
        newE = np.zeros((n, yrsOfAnalysis))
        for k in train_set:
            tmp_high_infected = normalized_high_infected_list[k]
            tmp_high_latent = normalized_high_latent_list[k]
            tmp_high_population = normalized_high_population_list[k]
            tmp_birth_rate = sfw.zeros((n, yrsOfAnalysis))
            tmp_newI = np.zeros((n, yrsOfAnalysis))
            tmp_newE = np.zeros((n, yrsOfAnalysis))
            for i in range(yrsOfAnalysis-1):
                tmp_birth_rate[0, i] = tmp_high_population[0,i+1]
                tmp_newI[0,i] = tmp_high_infected[0,i+1]
                tmp_newE[0,i] = tmp_high_latent[0,i+1]
            birth_rate += tmp_birth_rate
            newI += tmp_newI
            newE += tmp_newE
            #new_x, new_y = find_beta_constraints(yrsOfAnalysis, tmp_high_infected, tmp_high_population, tmp_newI, mu, nu_list[k], death_rate, tmp_birth_rate) # remember to modify "normalized..."
            new_x_infected, new_y_infected, new_x_latent, new_y_latent = util.separated_find_beta_constraints(yrsOfAnalysis, tmp_high_infected, tmp_high_population, tmp_high_latent, tmp_newI, tmp_newE, mu, nu_list[k], death_rate, tmp_birth_rate) # remember to modify "normalized..."
    
            x_train_infected_list.append(new_x_infected)
            y_train_infected_list.append(new_y_infected)
            x_train_latent_list.append(new_x_latent)
            y_train_latent_list.append(new_y_latent)

        birth_rate /= float(train_number)
        newI /= float(train_number)
        newE /= float(train_number)

        # ============================= training =============================
        print "======================== training ============================="
        x_train_infected = np.concatenate(x_train_infected_list, axis=1)
        y_train_infected = np.concatenate(y_train_infected_list, axis=1)
        x_train_latent = np.concatenate(x_train_latent_list, axis=1)
        y_train_latent = np.concatenate(y_train_latent_list, axis=1)
    
        x_train = np.concatenate([x_train_infected, x_train_latent], axis=1)
        y_train = np.concatenate([y_train_infected, y_train_latent] , axis=1)


        weight = np.concatenate( (np.ones((n, train_number * (yrsOfAnalysis - 1))) ,                  # infected part
                                  np.ones((n, train_number * (yrsOfAnalysis - 1)))/100     ), axis=1) # latent part


        print "====================== running linear regression ======================="
        beta, alpha = util.quadratic_solver(x_train, y_train, weight, yrsOfAnalysis)

        # =============== optimization of low fidelity model =================

        L = nu_list[-1] * 0.8
        U = np.ones(110) * 0.3
        K = effort - np.sum(L)
        num_iter = 100
        c = np.bmat([[np.ones((n, 1))], [np.zeros((n+1, 1))]])

        nu = L
        for i in range(num_iter):
            print i

            x0 = np.transpose(np.bmat([np.transpose(I_infected[:,0]), np.zeros((n)), np.ones((1))]))

#            grad = gradient(n, T, G, S[j], x0, c, newE, newI, nu+L, mu, d, alpha_fast, alpha_slow, beta[j], N[j])
            #Swap this line in for Tth step only, vs sum of 1...T 
#            grad = grad + gradient(n, T, G, S[j], x0, c, newE, newI, nu+L, mu, d, alpha_fast, alpha_slow, beta[j], N[j])

            #grad = sfw.gradient_sum(n, yrsOfAnalysis, G, S_safe, x0, c, newE, newI, nu+L, mu, death_rate, alpha_fast, alpha_slow, beta, N_population)

            derI, derE, derS = derivative.getDerivativeByNu(n, yrsOfAnalysis, I_infected, E_latent, S_safe, newI, birth_rate, death_rate, nu, mu, beta, alpha)
            grad = - np.sum(derI, axis=(1,2)) # TODO

            #grad = grad / len(beta)
            v = sfw.greedy(grad, U - L, np.zeros((n)), K)
            nu = nu + (1./num_iter)*v

            low_infected, low_population, low_latent = util.query_low(n, yrsOfAnalysis, G, S_safe, newE, newI, nu, mu, death_rate, alpha, beta, N_population, I_infected, E_latent, birth_rate)
            tmp_objective_value = np.sum(low_infected)
            print "objective value at round {0}: {1}".format(i, tmp_objective_value)

        # -------------------- varified the SFW algorithm --------------------
        for i in range(50):
            random_nu = util.get_random_nu(effort, n)

            low_infected, low_population, low_latent = util.query_low(n, yrsOfAnalysis, G, S_safe, newE, newI, random_nu, mu, death_rate, alpha, beta, N_population, I_infected, E_latent, birth_rate) # query_low(n, optTimeHorizon, G, S, newE, newI, nu, mu, d, alpha_fast, alpha_slow, beta, N, I, b):
            low_objective_value = np.sum(low_infected)
            print "random nu: {0}, objective value: {1}".format(i, low_objective_value)
            


        # ================ running both fidelity models to check =============
        # ----------------------- high fidelity model ------------------------

        #"""
        tmp_nu = np.resize(nu, (n,1))
        tmp_high_infected, tmp_high_latent, tmp_high_population = util.query_high(NumPpl, 25, np.concatenate([tmp_nu, tmp_nu], axis=1).tolist(), uptakePolicy.finerUptakeAgeBracs)

        tmp_high_infected = np.array(tmp_high_infected)
        tmp_high_population = np.array(tmp_high_population)
        tmp_high_latent = np.array(tmp_high_latent)

        normalized_high_infected = tmp_high_infected[:, range(130*12, 155*12, 12)]
        normalized_high_population = tmp_high_population[:, range(130*12, 155*12, 12)]
        normalized_high_latent = tmp_high_latent[:, range(130*12, 155*12, 12)]

        normalized_high_infected_list.append(normalized_high_infected)
        normalized_high_population_list.append(normalized_high_population)
        normalized_high_latent_list.append(normalized_high_latent)

        nu_list.append(nu)

        high_objective_value = np.sum(normalized_high_infected) 

        # ----------------------- low fidelity model -------------------------
        low_infected, low_population, low_latent = util.query_low(n, yrsOfAnalysis, G, S_safe, newE, newI, nu, mu, death_rate, alpha, beta, N_population, I_infected, E_latent, birth_rate) # query_low(n, optTimeHorizon, G, S, newE, newI, nu, mu, d, alpha_fast, alpha_slow, beta, N, I, b):
        low_objective_value = np.sum(low_infected)

        # ----------------------------- recording ----------------------------
        low_objective_value_list.append(low_objective_value)
        high_objective_value_list.append(high_objective_value)

        print "iteration: {0}, high objective value: {1}, low objective value: {2}".format(iteration, high_objective_value, low_objective_value)
        infected_variation = np.sum(np.abs(normalized_high_infected - low_infected))
        population_variation = np.sum(np.abs(normalized_high_population - low_population))
        print "total variation infected: {0}, population: {1}".format(infected_variation, population_variation)


        #"""

        # ============================= pickle save ==============================
        pickle.dump((normalized_high_infected_list, normalized_high_latent_list, normalized_high_population_list, nu_list), open(dir_path + "data/high_fidelity_{0}_{1}_{2}.data".format(date, file_index, NumPpl), "wb"))
