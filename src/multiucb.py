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
import sklearn
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel, Matern
from sklearn import linear_model
from sklearn.model_selection import cross_val_predict
import scipy
import copy

import params
import util
import customized_gpr

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

def decomposed_mean_std(gpr_list, x):
    mean = 0
    variance = 0
    for k in range(len(gpr_list)):
        tmp_mean, tmp_std = gpr_list[k].predict(x.reshape(1,-1), return_std=True)
        mean += tmp_mean
        variance += tmp_std**2
    std = np.sqrt(variance)

    return mean, std

def decomposed_GPUCB_ObjectiveValue(gpr_list, x, t, n, delta=0.1, b=1, a=1, r=1):
    #bt = 2 * np.log(t**2 * 2 * np.pi**2 / (3 * delta)) + 2 * n * np.log(t**2 * n * b * r * np.sqrt(np.log(4 * n * a / delta)))
    bt = 2 * np.log(t**2 * 2 * np.pi**2 / (3 * delta)) + 2 * n * np.log(t**2 * n * b * r * np.sqrt(np.log(4 * n * a / delta))) / 400
    mean = 0
    variance = 0
    for k in range(len(gpr_list)):
        tmp_mean, tmp_std = gpr_list[k].predict(x.reshape(1,-1), return_std=True)
        mean += tmp_mean
        variance += tmp_std**2
    std = np.sqrt(variance)

    return mean - np.sqrt(bt) * std

def bilinear_constraints(yrsOfAnalysis, high_infected, high_population, high_latent, mu, nu, death_rate, discount_factor=0.5): # normalized high_infected and high_population
    MAX_AGE = 110

    nu = nu * discount_factor

    high_healthy = high_population - high_infected - high_latent
    mu = np.array(mu.tolist())

    # ====================== infected regression ======================
    new_x_infected_list = np.zeros((yrsOfAnalysis-1, MAX_AGE, 1, MAX_AGE + 2)) # plus alpha_i and constant
    new_y_infected_list = np.zeros((yrsOfAnalysis-1, MAX_AGE, 1))
    new_x_latent_list = np.zeros((yrsOfAnalysis-1, MAX_AGE, 1, MAX_AGE + 2))
    new_y_latent_list = np.zeros((yrsOfAnalysis-1, MAX_AGE, 1))
    new_x_healthy_list = np.zeros((yrsOfAnalysis-1, MAX_AGE, 1, MAX_AGE + 2))
    new_y_healthy_list = np.zeros((yrsOfAnalysis-1, MAX_AGE, 1))

    for t in range(yrsOfAnalysis - 1):
        for i in range(MAX_AGE-1):
            # ------------------- y target --------------------
            new_y_latent_list[t][i][0] = high_latent[i+1][t+1] - high_latent[i][t] * (1-mu[i][t]) # - newE[i+1][t]
            new_y_infected_list[t][i][0] = high_infected[i+1][t+1] - high_infected[i][t] * (1 - death_rate[i]) * (1 - nu[i]) # - newI[i+1][t]
            new_y_healthy_list[t][i][0] = high_healthy[i+1][t+1] - high_infected[i][t] * (1 - death_rate[i]) * nu[i]

            # ------------- new exposed latent ----------------
            for k in range(MAX_AGE):
                if high_population[k][t] > 0:
                    new_exposed = high_healthy[i][t] * (1 - mu[i][t]) * (float(high_infected[k][t]) / high_population[k][t])
                    new_x_latent_list[t][i][0][k] = new_exposed
                    new_x_healthy_list[t][i][0][k] = - new_exposed
            new_x_latent_list[t][i][0][MAX_AGE + 0] = - high_latent[i][t] * (1 - mu[i][t]) # alpha coefficient
            new_x_infected_list[t][i][0][MAX_AGE + 0] = high_latent[i][t] * (1 - mu[i][t]) # alpha coefficient

            # ---------------- constant part ------------------
            new_x_latent_list[t][i][0][MAX_AGE + 1] = 1
            new_x_infected_list[t][i][0][MAX_AGE + 1] = 1
            new_x_healthy_list[t][i][0][MAX_AGE + 1] = 1

    return new_x_infected_list, new_x_latent_list, new_x_healthy_list, new_y_infected_list, new_y_latent_list, new_y_healthy_list

def qp_solver(x_train, y_train, age_group_index, beta_nearby_difference=0.1):
    training_number, variable_number = x_train.shape
    assert(variable_number == MAX_AGE + 2)
    
    model = cplex.Cplex()
    model.set_log_stream(None)
    model.set_error_stream(None)
    model.set_warning_stream(None)
    model.set_results_stream(None)

    model.objective.set_sense(model.objective.sense.minimize)

    x_valid = np.matrix(x_train)
    y_valid = np.matrix(y_train.reshape(-1,1))
    x_valid_sum = np.sum(np.abs(np.array(x_valid)), axis=0)
    x_valid_nonzero = np.logical_or((x_valid_sum > 0.01), (x_valid_sum < -0.01))

    quadratic_matrix = np.array(np.transpose(x_valid) * x_valid) * 2

    linear_coefficients = np.array(-2 * np.transpose(x_valid) * y_valid)[:,0]
    constant = float(np.transpose(y_valid) * y_valid)

    ub = np.ones(variable_number)
    lb = np.zeros(variable_number)
    # -------------- coefficients ------------------
    model.variables.add(obj=linear_coefficients, ub=ub, lb=lb, names=["beta_{0}_{1}".format(i, j) for j in range(MAX_AGE + 1)] + ["alpha_{0}".format(i)])

    qmat = []
    for j in range(variable_number):
        qmat.append([range(variable_number), list(quadratic_matrix[j]) ])

    model.objective.set_quadratic(qmat)

    ## ============== transition matrix beta ============= #
    for j in range(MAX_AGE):
        if j != age_group_index:
            model.linear_constraints.add(lin_expr=[[[age_group_index, j], [1, -1]]], senses=["G"], rhs=[0])
        if j+1 < MAX_AGE:
            model.linear_constraints.add(lin_expr=[[[j, j+1], [1, -1]]], senses=["L"], rhs=[beta_nearby_difference])
            model.linear_constraints.add(lin_expr=[[[j+1, j], [1, -1]]], senses=["L"], rhs=[beta_nearby_difference])

    model.solve()
    objective_value = model.solution.get_objective_value() + constant
    #print("age: {0}, obj: {1}".format(age_group_index, objective_value))
    solutions = model.solution.get_values() * x_valid_nonzero

    return solutions


def linear_predict(x_train, y_train, solutions):
    prediction = np.dot(x_train, solutions)
    error = y_train - prediction
    error_percentage = (error*100 / y_train) 

    print("error: {0}".format(error))
    print("error percentage: {0}".format(error_percentage))

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
                #effort = np.random.randint(10,30)
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
        #normalized_high_infected_list, normalized_high_latent_list, normalized_high_population_list, nu_list = pickle.load(open(dir_path + "data/high_fidelity_{0}_{1}_{2}.data".format(date, file_index, str(NumPpl)), "rb"), encoding='latin1')
        normalized_high_infected_list, normalized_high_latent_list, normalized_high_population_list, nu_list = pickle.load(open(dir_path + "data/high_fidelity_{0}_{1}_{2}.data".format(date, file_index, str(NumPpl)), "rb"))


    samples_number = len(normalized_high_infected_list)
    for i in range(samples_number):
        tmp_total_infected_population = np.sum(normalized_high_infected_list[i])
        objective_value_list.append(tmp_total_infected_population)


    # ===================== extracting the constraints ====================
    x_train_infected_list = []
    y_train_infected_list = []
    x_train_latent_list = []
    y_train_latent_list = []
    x_train_healthy_list = []
    y_train_healthy_list = []

    #samples_number = 10

    for i in range(samples_number):
        print("extracting constraint {0}...".format(i))
        new_x_infected, new_x_latent, new_x_healthy, new_y_infected, new_y_latent, new_y_healthy = bilinear_constraints(yrsOfAnalysis, normalized_high_infected_list[i], normalized_high_population_list[i], normalized_high_latent_list[i], mu, nu_list[i], death_rate)

        x_train_infected_list.append(new_x_infected)
        y_train_infected_list.append(new_y_infected)
        x_train_latent_list.append(new_x_latent)
        y_train_latent_list.append(new_y_latent)
        x_train_healthy_list.append(new_x_healthy)
        y_train_healthy_list.append(new_y_healthy)

    x_train_infected = np.concatenate(x_train_infected_list, axis=2)
    x_train_latent   = np.concatenate(x_train_latent_list, axis=2)
    x_train_healthy  = np.concatenate(x_train_healthy_list, axis=2)
    y_train_infected = np.concatenate(y_train_infected_list, axis=2)
    y_train_latent   = np.concatenate(y_train_latent_list, axis=2)
    y_train_healthy  = np.concatenate(y_train_healthy_list, axis=2)

    # =================== training linear regression =====================
    """
    for t in range(yrsOfAnalysis - 1):
        for i in range(MAX_AGE - 1):
            tmp_x_train_infected = x_train_infected[t][i] # / (np.mean(y_train_infected[t][i]) + 1)
            tmp_y_train_infected = y_train_infected[t][i] # / (np.mean(y_train_infected[t][i]) + 1)
            tmp_x_train_latent   = x_train_latent[t][i]  / 10 # / (np.mean(y_train_latent[t][i])   + 1) 
            tmp_y_train_latent   = y_train_latent[t][i]  / 10 # / (np.mean(y_train_latent[t][i])   + 1)
            tmp_x_train_healthy  = x_train_healthy[t][i] / 1000 # / (np.mean(y_train_healthy[t][i])  + 1)
            tmp_y_train_healthy  = y_train_healthy[t][i] / 1000 # / (np.mean(y_train_healthy[t][i])  + 1)


            tmp_x_train = np.concatenate([tmp_x_train_infected, tmp_x_train_latent, tmp_x_train_healthy], axis=0)
            tmp_y_train = np.concatenate([tmp_y_train_infected, tmp_y_train_latent, tmp_y_train_healthy], axis=0)
            lr = linear_model.LinearRegression()

            # ------------------ linear regression model -------------------
            #lr.fit(tmp_x_train, tmp_y_train)
            #y_predict = cross_val_predict(lr, tmp_x_train, tmp_y_train, cv=10)
            #R_score = lr.score(tmp_x_train, tmp_y_train)

            # -------------------- quadratic programming -------------------
            # this model provides a more reasonable explanation but more limited
            solutions = qp_solver(tmp_x_train, tmp_y_train, i)
            y_predict = np.dot(tmp_x_train, solutions)

            error = y_predict - tmp_y_train
            score = np.std(error)
            mean = np.mean(error)
            y_mean = np.mean(tmp_y_train)
            y_std = np.std(tmp_y_train)
            print("linear regression at t: {0}, i: {1}, std score: {2}, mean: {3}, y mean: {4}, y std: {5}".format(t, i, score, mean, y_mean, y_std))
    #"""


    # ======================== GP UCB algorithm ============================
    total_iterations = 100
    for iteration_index in range(len(nu_list)+1, total_iterations+len(nu_list)+1):
        # ------------------------- GP regression --------------------------
        objective_value_std = np.std(objective_value_list)
        objective_value_mean = np.mean(objective_value_list)
        X_train = np.array(nu_list)
        #y_train = (np.array(objective_value_list) - objective_value_mean) / objective_value_std


        # ------------------------ sklearn gp regression -------------------
        gpr_list = []
        mean_list = []
        std_list = []
        #tmp_X_train = X_train
        #universal_gpr = customized_gpr.GaussianProcessRegressor()
        #universal_gpr.fit_X(tmp_X_train)
        for t in range(yrsOfAnalysis):
        #for i in range(MAX_AGE):
            #print("gaussian process regression at time {0}, age {1}".format(t, i))
            print("gaussian process regression at time {0}".format(t))
            #training_size = int(np.ceil(len(X_train) * 0.8))
            #training_indices = np.random.choice(len(nu_list), training_size, replace=False)
            tmp_X_train = X_train
            tmp_y_train = np.sum(np.array(normalized_high_infected_list)[:,:,t], axis=1)
            #tmp_y_train = np.array(normalized_high_infected_list)[:,i,t]
            #tmp_y_mean = np.mean(tmp_y_train)
            #tmp_y_std = np.std(tmp_y_train)
            #if tmp_y_std != 0:
            #    tmp_y_train = (tmp_y_train - tmp_y_mean) / tmp_y_std

            #"""
            #kernel = ConstantKernel(1, constant_value_bounds=(1e-4,1)) * RBF(0.1, (1e-2, 1)) + WhiteKernel(noise_level=1e-4, noise_level_bounds=(1e-7,1e-3))
            kernel = ConstantKernel() * RBF()
            #kernel = ConstantKernel(1, constant_value_bounds=(1e-2,2e0)) * Matern(length_scale=0.1, nu=3.0/2) + WhiteKernel(noise_level=0.5, noise_level_bounds=(1e-3,1))
            gpr = GaussianProcessRegressor(kernel=kernel, alpha=1e-2, normalize_y=True)
            #gpr = GaussianProcessRegressor() # no kernel
            gpr.fit(tmp_X_train, tmp_y_train)
            #"""


            """
            gpr = copy.copy(universal_gpr)
            gpr.fit_y(tmp_y_train)
            #"""

            gpr_list.append(gpr)
            #mean_list.append(tmp_y_mean)
            #std_list.append(tmp_y_std)
            
        fn = lambda x: decomposed_GPUCB_ObjectiveValue(gpr_list, x, iteration_index, n)

        # ------------------------- gpflow regression ----------------------
        #kernel = gpflow.kernels.RBF(1) * gpflow.kernels.Constant(100)
        #gpr = gpflow.models.gpr.GPR(X_train, y_train, kern=kernel)
        #fn = lambda x: GPUCB_ObjectiveValue(gpr, x, t, n)

        # ---------------------------- minimization ------------------------
        print("minimization...")

        cons = ({"type": "eq", "fun": lambda x: np.sum(x) - budget_constraint})
        bds = [(0, 1) for i in range(n)]
        #initial_nu = np.array(nu_list[0])
        initial_nu = np.zeros(n)

        res = scipy.optimize.minimize(fn, initial_nu.reshape(n,1), bounds=bds, constraints=cons)

        new_nu = res.x
        objective_value = res.fun
        #mean, std = gpr.predict(new_nu.reshape(1,-1), return_std=True)
        mean, std = decomposed_mean_std(gpr_list, new_nu)
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

        print("iteration: {0}, mean: {1}, std: {2}, mean - std: {3}, true outcome: {4}".format(iteration_index, mean, std, mean-std, total_infected_population))

        pickle.dump((normalized_high_infected_list, normalized_high_latent_list, normalized_high_population_list, nu_list), open(dir_path + "data/high_fidelity_{0}_{1}_{2}.data".format(date, file_index, str(NumPpl)), "wb"), protocol=2)

