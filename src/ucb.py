import matlab.engine
import os
import uptakePolicy
import numpy as np
import itertools
import pickle
import sys
import argparse
#import cplex
import derivative
import sklearn
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel, Matern
import scipy

import params
import util

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

# ================ functions ==================

def GPUCB_ObjectiveValue(gpr, x, t, n, delta=0.1, b=1, a=1, r=1):
    bt = 2 * np.log(t**2 * 2 * np.pi**2 / (3 * delta)) + 2 * n * np.log(t**2 * n * b * r * np.sqrt(np.log(4 * n * a / delta)))
    #bt = 2 * np.log(t**2 * 2 * np.pi**2 / (3 * delta)) + 2 * n * np.log(t**2 * n * b * r * np.sqrt(np.log(4 * n * a / delta))) / 400
    mean, std = gpr.predict(x.reshape(1,-1), return_std=True)

    return mean - np.sqrt(bt) * std

def customized_RBF(x, y, length_scale=1.0):
    return np.exp(- 0.5 * (np.dot(x - y, x - y)) / length_scale**2)

def customized_RBF_derivative(x, y, length_scale=1.0): # d k(x,y)/ d x
    return - np.exp(- 0.5 * (np.dot(x - y, x - y)) / length_scale**2) * (x - y) / (length_scale ** 2)



if __name__ == "__main__":
    restart = False
    restart_training = True
    rerun = False   # high fidelity model
    enable_plot = False
    samples_number = 1
    budget_constraint = 10

    disease = 'tb'

    seed = 1234
    np.random.seed(seed)

    #NumPpl = 20000
    yrsOfAnalysis = 25
    n = 110

    parser = argparse.ArgumentParser(description="filename parser")
    parser.add_argument("-d", "--date", help="Input date")
    parser.add_argument("-m", "--method", help="Input one of the methods: initial, optimal, full, peak40, peak55")
    parser.add_argument("-n", "--number", help="Input the number of initial population")
    parser.add_argument("-v", "--variance", default=0, help="Input the change to the rerun nu")
    parser.add_argument("-t", "--testingSample", default=0, help="Input the number of testing sample")

    args = parser.parse_args()
    print("args: {0}".format(args))

    date = args.date
    file_index = args.method
    NumPpl = int(args.number)
    rerun_variance = float(args.variance)
    testing_sample = int(args.testingSample)

    print("date: {0}, number of population {1}".format(date, str(NumPpl)))

    f, p, _, _, _ = sfw.load_data(disease, mat_path=Bryan_path + "SIS_forBryan_TB.mat")

    # =================== parameter setting ===================
    # ----------- precomputed low fidelity setting ------------

    I_infected = np.zeros((n, yrsOfAnalysis))
    N_population = np.zeros((n, yrsOfAnalysis))
    E_latent = np.zeros((n, yrsOfAnalysis))
    S_safe = np.zeros((n, yrsOfAnalysis))

    beta = np.zeros((n,n))

    # -------------------- default parameters ------------------
    newE, newI, alpha_fast, alpha_slow, G, mu, death_rate = util.get_parameters(n, yrsOfAnalysis)

    # ==================== precomputed data ====================
    # =================== high fidelity model ==================

    uptakeUrbanKnowledge = np.array(uptakePolicy.finerUptakeUrbanKnowledge)

    nu_list = []
    objective_value_list = []
    normalized_high_infected_list = []
    normalized_high_population_list = []
    normalized_high_latent_list = []

    if restart:
        for i in range(samples_number):
            if file_index == "random":
                effort = 10
                nu = util.get_random_nu(effort, n)
            nu_list.append(nu)
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

        pickle.dump((normalized_high_infected_list, normalized_high_latent_list, normalized_high_population_list, nu_list), open(dir_path + "data/high_fidelity_{0}_{1}_{2}.data".format(date, file_index, str(NumPpl)), "wb"), protocol=2)
    else:
        normalized_high_infected_list, normalized_high_latent_list, normalized_high_population_list, nu_list = pickle.load(open(dir_path + "data/high_fidelity_{0}_{1}_{2}.data".format(date, file_index, str(NumPpl)), "rb"), encoding='latin1')


    samples_number = len(normalized_high_infected_list)
    for i in range(samples_number):
        tmp_total_infected_population = np.sum(normalized_high_infected_list[i])
        objective_value_list.append(tmp_total_infected_population)


    # ======================== GP UCB algorithm ============================
    total_iterations = 100
    for t in range(len(nu_list)+1, total_iterations+len(nu_list)+1):
        # ------------------------- GP regression --------------------------
        objective_value_std = np.std(objective_value_list)
        objective_value_mean = np.mean(objective_value_list)
        X_train = np.array(nu_list)
        #y_train = (np.array(objective_value_list) - objective_value_mean) / objective_value_std
        y_train = np.array(objective_value_list)


        # ------------------------ sklearn gp regression -------------------
        #kernel = ConstantKernel(1, constant_value_bounds=(1e-4,1)) * RBF(0.1, (1e-2, 1)) + WhiteKernel(noise_level=1e-4, noise_level_bounds=(1e-7,1e-3))
        kernel = ConstantKernel() * RBF()
        #kernel = ConstantKernel(1, constant_value_bounds=(1e-2,2e0)) * Matern(length_scale=0.1, nu=3.0/2) + WhiteKernel(noise_level=0.5, noise_level_bounds=(1e-3,1))
        gpr = GaussianProcessRegressor(alpha=0.25, normalize_y=True ,kernel=kernel)
        #gpr = GaussianProcessRegressor() # no kernel
        gpr.fit(X_train, y_train)
        fn = lambda x: GPUCB_ObjectiveValue(gpr, x, t, n)

        # ------------------------- gpflow regression ----------------------
        #kernel = gpflow.kernels.RBF(1) * gpflow.kernels.Constant(100)
        #gpr = gpflow.models.gpr.GPR(X_train, y_train, kern=kernel)
        #fn = lambda x: GPUCB_ObjectiveValue(gpr, x, t, n)

        # ---------------------------- minimization ------------------------

        cons = ({"type": "eq", "fun": lambda x: np.sum(x) - budget_constraint})
        bds = [(0, 1) for i in range(n)]
        #initial_nu = np.array(nu_list[0])
        initial_nu = np.zeros(n)

        res = scipy.optimize.minimize(fn, initial_nu.reshape(n,1), bounds=bds, constraints=cons)

        new_nu = res.x
        objective_value = res.fun * objective_value_std + objective_value_mean
        mean, std = gpr.predict(new_nu.reshape(1,-1), return_std=True)
        #std = std * objective_value_std
        #mean = mean * objective_value_std + objective_value_mean

        # ----------------------- run high fidelity model ------------------
        tmp_nu = np.resize(new_nu, (n,1))
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
 
        nu_list.append(new_nu)
        total_infected_population = np.sum(normalized_high_infected)
        objective_value_list.append(total_infected_population)

        print("iteration: {0}, mean: {1}, std: {2}, mean - std: {3}, true outcome: {4}".format(t, mean, std, mean-std, total_infected_population))

        pickle.dump((normalized_high_infected_list, normalized_high_latent_list, normalized_high_population_list, nu_list), open(dir_path + "data/high_fidelity_{0}_{1}_{2}.data".format(date, file_index, str(NumPpl)), "wb"), protocol=2)

