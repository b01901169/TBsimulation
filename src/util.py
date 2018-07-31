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

# ================ functions ==================

# ------------- copied from Bryan -------------
def single_infected_sim(n, T, G, S_safe, newE, newI, nu, mu, death_rate, alpha, single_beta, N_population, I_infected, E_latent, birth_rate):
    import numpy as np
    I_infected = np.matrix(I_infected)
    N_population = np.matrix(N_population)
    E_latent = np.matrix(E_latent)
    S_safe = np.matrix(S_safe)

    I = sfw.zeros((n, T))
    I[:, 0] = I_infected[:,0]
    N = sfw.zeros((n, T))
    N[:, 0] = N_population[:, 0]
    S = sfw.zeros((n, T))
    S[:, 0] = S_safe[:, 0]
    E = sfw.zeros((n, T))
    E[:, 0] = E_latent[:, 0]
    beta = single_beta
    for t in range(1, T):
        infected_proportion = np.ones(n)/ np.resize(np.array(N[:,t-1]), n)
        infected_proportion[infected_proportion == np.inf] = 0
        newExposed = sfw.diag(np.resize(S[:,t-1], (n))) * sfw.diag(1 - np.resize(mu[:,t-1], n)) * beta * sfw.diag(infected_proportion) * I[:,t-1]
        #print newExposed
        #print N[:,t-1]
        E[:, t] = G*(sfw.diag(1-alpha[:,t-1] )*sfw.diag(1-mu[:,t-1])*E[:,t-1] + newExposed) + np.transpose(np.matrix(newE[:, t-1]))
        I[:, t] = G*(sfw.diag(1-death_rate) * sfw.diag(1-nu) * I[:,t-1] + sfw.diag(alpha[:,t-1] ) * sfw.diag(1-mu[:,t-1]) * E[:,t-1]) + np.transpose(np.matrix(newI[:,t-1]))
        S[:, t] = np.matrix(np.resize(birth_rate[:,t-1], (n,1))) + G*(sfw.diag(1-np.resize(mu[:,t-1], n))*S[:,t-1] + sfw.diag(nu)*sfw.diag(1-death_rate)*I[:,t-1] - newExposed)
        N[:, t] = S[:, t] + I[:, t] + E[:, t]
    return I, N, E

def query_low(n, optTimeHorizon, G, S, newE, newI, nu, mu, d, alpha, beta, N, I, E, b):
    disease = 'tb'
    #f, p, beta, n, G = sfw.load_data(disease, mat_path=Bryan_path + "SIS_forBryan_TB.mat")
    c = np.bmat([[sfw.ones((n, 1))], [sfw.zeros((n+1, 1))]])

    #if disease == 'gon':
    #    L = p.nu_sampled[1][:, p.T-1]
    #elif disease == 'tb':
    #    L = p.nu[:,p.T-1]
    #else:
    #    raise Exception('bad disease name')
    #Ks = [Ks[1], Ks[7]]
#    params = itertools.product(range(1, f.T+1), np.linspace(0.05, 0.4, 10))

    #print "p.T: {0}".format(p.T)
    #print "f.N: {0}".format(f.N)
    #print "L: {0}".format(L)
    #nu = L

    #infected_nu = sfw.infected_sim(n, optTimeHorizon, G, f.S, c, f.newE, f.newI, nu, f.mu, f.d, f.alpha_fast, f.alpha_slow, beta, f.N, f.I, f.b)
    #infected_nu, N_sim = sfw.single_infected_sim(n, optTimeHorizon, G, f.S[0], c, f.newE, f.newI, nu, f.mu, f.d, f.alpha_fast, f.alpha_slow, beta[0], f.N[0], f.I[0], f.b)

    #infected_nu, N_sim = sfw.single_infected_sim(n, optTimeHorizon, G, S, c, newE, newI, nu, mu, d, alpha_fast, alpha_slow, beta, N, I, b)
    I_sim, N_sim, E_sim = single_infected_sim(n, optTimeHorizon, G, S, newE, newI, nu, mu, d, alpha, beta, N, I, E, b)

    I_sim = np.array(I_sim)
    N_sim = np.array(N_sim)
    E_sim = np.array(E_sim)

    #print infected_nu, N_sim
    return I_sim, N_sim, E_sim


def query_high(NumPpl, yrsOfAnalysis, uptakeUrbanKnowledge, uptakeAgeBracs):
    eng_path = "/home/kai/Dropbox/USC/publication/TBsimulationCode/TBsimulationCode"
    eng = matlab.engine.start_matlab()
    eng.cd(eng_path)
    
    Group1, Group1_latent, age_population = eng.TBsimulation_jan23('.','p01', 130+yrsOfAnalysis, NumPpl, '-r0','loadBurnIn', 2018,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3,'base', 'customized',
                                    matlab.double(uptakeUrbanKnowledge), matlab.int8(uptakeAgeBracs), nargout=3)

    return Group1, Group1_latent, age_population

def get_random_nu(total_resource, MAX_AGE):
    max_number = 10
    min_number = 2
    initial_random_nu = np.random.randint(min_number, max_number, MAX_AGE)
    initial_sum = sum(initial_random_nu)
    nu = total_resource * initial_random_nu / float(initial_sum)
    return nu

def find_beta_constraints(yrsOfAnalysis, high_infected, high_population, newI, mu, nu, death_rate, birth_rate): # normalized high_infected and high_population
    high_healthy = high_population - high_infected
    mu = np.array(mu.tolist())
    birth_rate = np.array(birth_rate.tolist())

    # ====================== infected regression ======================
    new_x1 = np.zeros(((yrsOfAnalysis-1)*(MAX_AGE-1), MAX_AGE*MAX_AGE))
    new_y1 = np.zeros((yrsOfAnalysis-1)*(MAX_AGE-1))

    new_x2 = np.zeros(((yrsOfAnalysis-1)*(MAX_AGE-1), MAX_AGE*MAX_AGE))
    new_y2 = np.zeros((yrsOfAnalysis-1)*(MAX_AGE-1))
    for t in range(0, yrsOfAnalysis-1):
        for i in range(MAX_AGE-1):
            data_index = t * (MAX_AGE-1) + i
            # -------- infected people ---------
            new_y1[data_index] = high_infected[i+1][t+1] - (high_infected[i][t]) * (1-nu[i]) * (1-death_rate[i]) - newI[i+1][t]
            # --------- healthy people ---------
            new_y2[data_index] = high_healthy[i+1][t+1] - birth_rate[i+1][t] - nu[i] * (1-death_rate[i]) * high_infected[i][t]

            # ----------- new exposed ---------- # haven't considered alpha_fast
            for k in range(MAX_AGE):
                beta_index = i*MAX_AGE + k
                if high_population[k][t] > 0:
                    new_exposed = high_healthy[i][t] * (1 - mu[i][t]) * (float(high_infected[k][t]) / high_population[k][t])
                    new_x1[data_index][beta_index] = new_exposed
                    new_x2[data_index][beta_index] = -new_exposed
    new_x = np.concatenate((new_x1, new_x2))
    new_y = np.concatenate((new_y1, new_y2))

    #return new_x1, new_y1
    return new_x, new_y

def separated_find_beta_constraints(yrsOfAnalysis, high_infected, high_population, high_latent, newI, newE, mu, nu, death_rate, birth_rate): # normalized high_infected and high_population
    MAX_AGE = 110

    high_healthy = high_population - high_infected - high_latent
    mu = np.array(mu.tolist())
    birth_rate = np.array(birth_rate.tolist())

    # ====================== infected regression ======================
    new_x_infected_list = np.zeros((MAX_AGE, yrsOfAnalysis-1, MAX_AGE + yrsOfAnalysis - 1)) # plus alpha_i
    new_y_infected_list = np.zeros((MAX_AGE, yrsOfAnalysis-1))
    new_x_latent_list = np.zeros((MAX_AGE, yrsOfAnalysis-1, MAX_AGE + yrsOfAnalysis - 1))
    new_y_latent_list = np.zeros((MAX_AGE, yrsOfAnalysis-1))

    for i in range(MAX_AGE-1):
        for t in range(0, yrsOfAnalysis-1):
            # -------- infected people ---------
            new_y_latent_list[i][t] = high_latent[i+1][t+1] - high_latent[i][t] * (1-mu[i][t]) - newE[i+1][t]
            new_y_infected_list[i][t] = high_infected[i+1][t+1] - high_infected[i][t] * (1 - death_rate[i]) * (1 - nu[i]) - newI[i+1][t]

            # ----------- new exposed ---------- # haven't considered alpha_fast
            for k in range(MAX_AGE):
                if high_population[k][t] > 0:
                    new_exposed = high_healthy[i][t] * (1 - mu[i][t]) * (float(high_infected[k][t]) / high_population[k][t])
                    new_x_latent_list[i][t][k] = new_exposed
            new_x_latent_list[i][t][MAX_AGE + t] = - high_latent[i][t] * (1 - mu[i][t]) # alpha coefficient
            new_x_infected_list[i][t][MAX_AGE + t] = high_latent[i][t] * (1 - mu[i][t]) # alpha coefficient

    #return_x_list = np.concatenate([new_x_latent_list, new_x_infected_list], axis=1)
    #return_y_list = np.concatenate([new_y_latent_list, new_y_infected_list], axis=1)

    return new_x_infected_list, new_y_infected_list, new_x_latent_list, new_y_latent_list
    #return return_x_list, return_y_list

def linear_regression(x_train, y_train):
    from sklearn import linear_model

    # Create linear regression object
    regr = linear_model.LinearRegression()

    regr.fit(x_train, y_train)
    coefficients = regr.coef_
    score = regr.score(x_train, y_train)

    print("Score: {0}".format(score))
    print("Coefficients: \n", coefficients)
    new_beta = np.resize(coefficients, (110, 110))
    return new_beta

def separated_linear_regression(x_train, y_train):
    from sklearn import linear_model

    # Initialization
    new_beta = np.zeros((MAX_AGE, MAX_AGE))

    # Create linear regression object
    for i in range(MAX_AGE):
        regr = linear_model.LinearRegression()
    
        regr.fit(x_train[i], y_train[i])
        coefficients = regr.coef_
        score = regr.score(x_train[i], y_train[i])
    
        print("{0}th Score: {1}".format(i, score))
        new_beta[i] = coefficients
    return new_beta

def quadratic_solver(x_train, y_train, weight, yrsOfAnalysis, beta_nearby_difference=0.1, alpha_nearby_difference=0):
    _, training_number, number_of_variables = x_train.shape
    assert(weight.shape[0] == MAX_AGE and weight.shape[1] == training_number)
    weight = np.sqrt(weight)

    new_beta = np.zeros((MAX_AGE, MAX_AGE))
    new_alpha = np.zeros((MAX_AGE, yrsOfAnalysis-1))

    for i in range(MAX_AGE-1):
        # ------------------ precompute coefficients ---------------
        weight_valid = np.matrix(np.diag(weight[i]))
        x_valid = np.matrix(x_train[i])
        y_valid = np.matrix(np.resize(y_train[i], (training_number, 1)))
        x_valid_sum = np.sum(np.abs(np.array(x_valid)), axis=0)
        x_valid_nonzero = np.logical_or((x_valid_sum > 0.01), (x_valid_sum < -0.01))
        #for j in range(MAX_AGE, MAX_AGE + yrsOfAnalysis - 1):
        #    x_valid_nonzero[j] = True    # alpha term

        x_weight = weight_valid * x_valid
        y_weight = weight_valid * y_valid

        quadratic_matrix = np.array(np.transpose(x_weight) * x_weight) * 2
 
        linear_coefficients = np.array(-2 * np.transpose(x_weight) * y_weight)[:,0]
        constant = float(np.transpose(y_weight) * y_weight)

        # -------------------- cplex setting ------------------------
        model = cplex.Cplex()
        model.set_log_stream(None)
        model.set_error_stream(None)
        model.set_warning_stream(None)
        model.set_results_stream(None)

        model.objective.set_sense(model.objective.sense.minimize)

        ub = np.ones(number_of_variables)
        lb = np.zeros(number_of_variables)
        # -------------- coefficients ------------------
        model.variables.add(obj=linear_coefficients, ub=ub, lb=lb, names=["beta_{0}_{1}".format(i, j) for j in range(number_of_variables-1)] + ["alpha_{0}".format(i)])
        # -------------- constant terms ----------------
        #model.variables.add()

        qmat = []
        for j in range(number_of_variables):
            qmat.append([range(number_of_variables), list(quadratic_matrix[j]) ])

        model.objective.set_quadratic(qmat)

        ## ============== transition matrix beta ============= #
        for j in range(MAX_AGE):
            if i != j:
                model.linear_constraints.add(lin_expr=[[[i,j], [1,-1]]], senses=["G"], rhs=[0])
            if j+1 < MAX_AGE:
                model.linear_constraints.add(lin_expr=[[[j, j+1], [1, -1]]], senses=["L"], rhs=[beta_nearby_difference])
                model.linear_constraints.add(lin_expr=[[[j+1, j], [1, -1]]], senses=["L"], rhs=[beta_nearby_difference])

        # ============== activation rate alpha ============== #
        for j in range(MAX_AGE, MAX_AGE + yrsOfAnalysis - 2):
            model.linear_constraints.add(lin_expr=[[[j, j+1], [1, -1]]], senses=["L"], rhs=[alpha_nearby_difference])
            model.linear_constraints.add(lin_expr=[[[j+1, j], [1, -1]]], senses=["L"], rhs=[alpha_nearby_difference])

        model.solve()
        objective_value = model.solution.get_objective_value() + constant
        print("age: {0}, obj: {1}".format(i, objective_value))
        solutions = model.solution.get_values() * x_valid_nonzero
        new_beta[i] = solutions[:MAX_AGE]
        new_alpha[i] = solutions[MAX_AGE:]

    return new_beta, new_alpha


def smart_linear_regression(x_train, y_train):
    MAX_AGE = 110
    x_delta = 2000
    y_delta = 0
    beta_list = np.zeros(MAX_AGE*MAX_AGE)
    large_matrix = np.zeros((MAX_AGE*MAX_AGE, MAX_AGE*MAX_AGE))
    for beta_index in range(MAX_AGE*MAX_AGE):
        valid_index = (x_train[:,beta_index] != 0)
        x_value = x_train[valid_index, beta_index]
        y_value = y_train[valid_index]
        x_sum = sum(x_value)
        y_sum = sum(y_value)
        if y_sum != 0:
            beta_value = (y_sum + y_delta) / float(x_sum + x_delta)
            beta_list[beta_index] = beta_value
    new_beta = np.resize(beta_list, (MAX_AGE,MAX_AGE))
    return new_beta

def get_parameters(n, yrsOfAnalysis):
    newE = np.zeros((n, yrsOfAnalysis))
    newI = np.zeros((n, yrsOfAnalysis))
    alpha_fast = np.ones(n)
    alpha_slow = np.ones(n)

    G = np.zeros((n,n))
    for i in range(n-1):
        G[i+1][i] = 1
    G = np.matrix(G)

    #death_rate = np.array([0.03003663] * n)    
    #mu = sfw.zeros((n, yrsOfAnalysis))
    mu, death_rate = params.getRateParameters(yrsOfAnalysis)
    return newE, newI, alpha_fast, alpha_slow, G, mu, death_rate


if __name__ == "__main__":
    restart = False
    restart_training = False
    rerun = False   # high fidelity model
    enable_plot = True
    samples_number = 20

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

    print("date: {0}, number of population {1}".format(date, NumPpl))

    print("loading ...")
    f, p, _, _, _ = sfw.load_data(disease, mat_path=Bryan_path + "SIS_forBryan_TB.mat")

    # ================== varied parameters ====================

    initial_nu = np.zeros(n)
    for i in range(n):
        modified_index = min(max(i-30,0), 30-1)
        initial_nu[i] = p.nu[modified_index,p.T-1]
    #initial_nu = np.array(p.nu[:,p.T-1]) # TODO

    # =================== parameter setting ===================
    # ----------- precomputed low fidelity setting ------------

    I_infected = np.zeros((n, yrsOfAnalysis))
    N_population = np.zeros((n, yrsOfAnalysis))
    E_latent = np.zeros((n, yrsOfAnalysis))
    S_safe = np.zeros((n, yrsOfAnalysis))

    beta = np.zeros((n,n))
    for brac_index1 in range(len(uptakePolicy.modTransMatAgeBrac)):
        for brac_index2 in range(len(uptakePolicy.modTransMatAgeBrac)):
            (start_age1, end_age1) = uptakePolicy.modTransMatAgeBrac[brac_index1]
            (start_age2, end_age2) = uptakePolicy.modTransMatAgeBrac[brac_index2]
            for i in range(start_age1, end_age1):
                for j in range(start_age2, end_age2):
                    beta[i][j] = uptakePolicy.modContactMat[brac_index1][brac_index2]

    initial_beta = beta
    initial_beta_sum = sum(sum(initial_beta))

    # -------------------- default parameters ------------------
    newE, newI, alpha_fast, alpha_slow, G, mu, death_rate = get_parameters(n, yrsOfAnalysis)

    #birth_rate = sfw.zeros((n, yrsOfAnalysis))
    initial_nu = np.zeros(n)
    for i in range(n):
        modified_index = min(max(i-30,0), 30-1)
        modified_birth_index = min(max(i-30,2), 30-1)
        #modified_index = i
        #mu[i] = f.mu[modified_index]
        initial_nu[i] = p.nu[modified_index,p.T-1]
    #for i in range(yrsOfAnalysis):
    #    birth_rate[0, i] = high_population[0,i]
        
    #d = np.array([f.d[0]] * n)
    #d = np.ones((n))/20


    # ==================== precomputed data ====================
    # =================== high fidelity model ==================
    
    uptakeUrbanKnowledge = np.array(uptakePolicy.finerUptakeUrbanKnowledge)

    nu_list = []
    normalized_high_infected_list = []
    normalized_high_population_list = []
    normalized_high_latent_list = []

    #nu = get_random_nu(10, n)

    if restart:
        for i in range(samples_number):
            if file_index == "random":
                effort = np.random.randint(10,30)
                nu = get_random_nu(effort, n)
                #nu[40:50] += 0.01
                #nu[30:40] -= 0.01
            nu_list.append(nu)
            tmp_nu = np.resize(nu, (n,1))
            tmp_high_infected, tmp_high_latent, tmp_high_population = query_high(NumPpl, 25, np.concatenate([tmp_nu, tmp_nu], axis=1).tolist(), uptakePolicy.finerUptakeAgeBracs)

            tmp_high_infected = np.array(tmp_high_infected)
            tmp_high_population = np.array(tmp_high_population)
            tmp_high_latent = np.array(tmp_high_latent)

            normalized_high_infected = tmp_high_infected[:, range(130*12, 155*12, 12)]
            normalized_high_population = tmp_high_population[:, range(130*12, 155*12, 12)]
            normalized_high_latent = tmp_high_latent[:, range(130*12, 155*12, 12)]

            normalized_high_infected_list.append(normalized_high_infected)
            normalized_high_population_list.append(normalized_high_population)
            normalized_high_latent_list.append(normalized_high_latent)

        pickle.dump((normalized_high_infected_list, normalized_high_latent_list, normalized_high_population_list, nu_list), open(dir_path + "data/high_fidelity_{0}_{1}_{2}.data".format(date, file_index, NumPpl), "wb"), protocol=2)
    else:
        normalized_high_infected_list, normalized_high_latent_list, normalized_high_population_list, nu_list = pickle.load(open(dir_path + "data/high_fidelity_{0}_{1}_{2}.data".format(date, file_index, NumPpl), "rb"), encoding="latin1")

    # =================== extract the constraints ==================
    x_train_infected_list = []
    y_train_infected_list = []
    x_train_latent_list = []
    y_train_latent_list = []

    samples_number = len(normalized_high_infected_list)
    #samples_number = 1

    birth_rate = sfw.zeros((n, yrsOfAnalysis))
    newI = np.zeros((n, yrsOfAnalysis))
    newE = np.zeros((n, yrsOfAnalysis))
    for k in range(samples_number):
        tmp_high_infected = normalized_high_infected_list[k]
        tmp_high_population = normalized_high_population_list[k]
        tmp_high_latent = normalized_high_latent_list[k]
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
        new_x_infected, new_y_infected, new_x_latent, new_y_latent = separated_find_beta_constraints(yrsOfAnalysis, tmp_high_infected, tmp_high_population, tmp_high_latent, tmp_newI, tmp_newE, mu, nu_list[k], death_rate, tmp_birth_rate) # remember to modify "normalized..."

        x_train_infected_list.append(new_x_infected)
        y_train_infected_list.append(new_y_infected)
        x_train_latent_list.append(new_x_latent)
        y_train_latent_list.append(new_y_latent)

    birth_rate /= float(samples_number)
    newI /= float(samples_number)
    newE /= float(samples_number)

    # ========== rerun the high fidelity model with new nu ===========
    print("=========== reruning the high fidelity model ===========")
    #testing_sample = 1
    if rerun:
        #nu = get_random_nu(10, n)
        #nu = nu_list[testing_sample]
        nu = np.ones(n)*10/70 # TODO
        nu[70:] = 0

        print("rerun variance: {0}".format(rerun_variance))
        nu[30:35] += rerun_variance
        nu[40:65] -= (rerun_variance / 5) # balance the deduction

        tmp_nu = np.resize(nu, (n,1))
        high_infected, high_latent, high_population = query_high(NumPpl, 25, np.concatenate([tmp_nu, tmp_nu], axis=1).tolist(), uptakePolicy.finerUptakeAgeBracs)

        high_infected = np.array(high_infected)
        high_population = np.array(high_population)
        high_latent = np.array(high_latent)

        high_infected = high_infected[:, range(130*12, 155*12, 12)]
        high_population = high_population[:, range(130*12, 155*12, 12)]
        high_latent = high_latent[:, range(130*12, 155*12, 12)]
    else:
        rerun_variance = 0
        nu = nu_list[testing_sample]
        high_infected = normalized_high_infected_list[testing_sample]
        high_population = normalized_high_population_list[testing_sample]
        high_latent = normalized_high_latent_list[testing_sample]

    high_percentage = np.zeros((n, yrsOfAnalysis))

    for i in range(n):
        for j in range(yrsOfAnalysis):
            if high_population[i][j] != 0:
                high_percentage[i][j] = float(high_infected[i][j]) / high_population[i][j]

    # ====================== linear regression =======================
    #x_train = np.concatenate(x_train_list)
    #y_train = np.concatenate(y_train_list)
    x_train_infected = np.concatenate(x_train_infected_list, axis=1)
    y_train_infected = np.concatenate(y_train_infected_list, axis=1)
    x_train_latent = np.concatenate(x_train_latent_list, axis=1)
    y_train_latent = np.concatenate(y_train_latent_list, axis=1)

    x_train = np.concatenate([x_train_infected, x_train_latent], axis=1)
    y_train = np.concatenate([y_train_infected, y_train_latent] , axis=1)

    weight = np.concatenate( (np.ones((n, samples_number * (yrsOfAnalysis - 1))) ,                 # infected part
                              np.ones((n, samples_number * (yrsOfAnalysis - 1)))/1000  ), axis=1)    # latent part
    #weight = np.ones((n, samples_number * (yrsOfAnalysis - 1)))
    #for i in range(n):
    #    for j in range(samples_number * (yrsOfAnalysis - 1)):
    #        weight[i][j] = float(high_population[i][j % (yrsOfAnalysis - 1)] ** 2)

    print("====================== running linear regression =======================")
    if restart_training:
        #new_beta = linear_regression(x_train, y_train)
        #new_beta = separated_linear_regression(x_train, y_train)
        new_beta, new_alpha = quadratic_solver(x_train, y_train, weight, yrsOfAnalysis)
        pickle.dump((new_beta, new_alpha), open(dir_path + "data/beta/beta_{0}_{1}_{2}.data".format(date, file_index, NumPpl), "wb"), protocol=2)
    else:
        new_beta, new_alpha = pickle.load( open(dir_path + "data/beta/beta_{0}_{1}_{2}.data".format(date, file_index, NumPpl), "rb"), encoding="latin1")

    #new_beta = np.zeros((MAX_AGE, MAX_AGE)) # TODO
    #new_alpha = np.zeros((MAX_AGE, yrsOfAnalysis-1)) # TODO
    #new_alpha = high_infected[:, 1:] / high_latent[:, :-1]
    #new_alpha[np.isnan(new_alpha)] = 0
    # ===================== plot setting ======================

    fig = plt.figure()
    beta_fig = plt.figure()

    # ------------------ plotting beta structure ----------------
    beta_ax = beta_fig.add_subplot(1,2,1)
    beta_ax.imshow(new_beta, cmap="Reds", interpolation="nearest", origin="lower")
    alpha_ax = beta_fig.add_subplot(1,2,2)
    alpha_ax.imshow(new_alpha, cmap="Reds", interpolation="nearest", origin="lower")

    # # ------------------- ploting 3D surface ------------------

    #"""
    #x_length, y_length = high_infected.shape
    x_start_year = 0 # e.g. 30 years old started
    x_end_year = 100 # e.g. 60 years old ended
    y_year_length = yrsOfAnalysis
    start_year = 1995
    x_high = [i for i in range(x_start_year, x_end_year) for j in range(y_year_length)]
    y_high = [start_year + j for i in range(x_start_year, x_end_year) for j in range(y_year_length)]
    data_length = (x_end_year - x_start_year) * (y_year_length)
    z_high_infected = [high_infected[x_high[i]][y_high[i] - start_year] for i in range(data_length)]
    z_high_population = np.array([high_population[x_high[i]][y_high[i] - start_year] for i in range(data_length)])
    z_high_latent = np.array([high_latent[x_high[i]][y_high[i] - start_year] for i in range(data_length)])
    z_high = np.array([high_percentage[x_high[i]][y_high[i] - start_year] for i in range(data_length)])

    ax1 = fig.add_subplot(2, 3, 1, projection='3d')
    ax1.set_xlabel("age")
    ax1.set_ylabel("year")
    ax1.set_zlabel("# latent")
     
    ax2 = fig.add_subplot(2, 3, 2, projection='3d')
    ax2.set_xlabel("age")
    ax2.set_ylabel("year")
    ax2.set_zlabel("# infected")
     
    ax3 = fig.add_subplot(2, 3, 3, projection='3d')
    ax3.set_xlabel("age")
    ax3.set_ylabel("year")
    ax3.set_zlabel("# population")

    ax1.plot_trisurf(x_high, y_high, z_high_latent, cmap=cm.coolwarm, linewidth=0.2, antialiased=True)
    ax2.plot_trisurf(x_high, y_high, z_high_infected, cmap=cm.coolwarm, linewidth=0.2, antialiased=True)
    ax3.plot_trisurf(x_high, y_high, z_high_population, cmap=cm.coolwarm, linewidth=0.2, antialiased=True)
    #"""

    # =================== parameter setting ===================
    # ----------- precomputed low fidelity setting ------------

    I_infected[:, 0] = high_infected[:,0]
    N_population[:, 0] = high_population[:,0]
    E_latent[:, 0] = high_latent[:,0]
    S_safe[:, 0] = N_population[:, 0] - I_infected[:, 0] - E_latent[:,0]

    birth_rate = np.zeros((n, yrsOfAnalysis))
    newI = np.zeros((n, yrsOfAnalysis))
    newE = np.zeros((n, yrsOfAnalysis))
    for i in range(yrsOfAnalysis-1):
        birth_rate[0, i] = high_population[0,i+1]
        newI[0,i] = high_infected[0,i+1]
        newE[0,i] = high_latent[0,i+1]

    beta = new_beta
    alpha = new_alpha

    # =================== low fidelity model ==================
    print("============================ low fidelity =============================")
    low_infected, low_population, low_latent = query_low(n, yrsOfAnalysis, G, S_safe, newE, newI, nu, mu, death_rate, alpha, beta, N_population, I_infected, E_latent, birth_rate) # query_low(n, optTimeHorizon, G, S, newE, newI, nu, mu, d, alpha_fast, alpha_slow, beta, N, I, b):

    pickle.dump((low_infected, low_population), open(dir_path + "data/low_fidelity_{0}_{1}.data".format(date, file_index), "wb"))

    # ------------------- ploting 3D surface ------------------
    #low_infected = np.array(low_infected)
    #low_population = np.array(low_population)
    low_percentage = np.zeros((n, yrsOfAnalysis))
    for i in range(n):
        for j in range(yrsOfAnalysis):
            if low_population[i][j] != 0:
                low_percentage[i][j] = float(low_infected[i][j]) / low_population[i][j]

    x_low = [i for i in range(x_start_year, x_end_year) for j in range(y_year_length)]
    y_low = [j + start_year for i in range(x_start_year, x_end_year) for j in range(y_year_length)]
    data_length = (x_end_year - x_start_year) * y_year_length
    low_infected = np.asarray(low_infected)
    z_low_infected = [low_infected[x_low[i]][y_low[i] - start_year] for i in range(data_length)]
    z_low_population = np.array([low_population[x_low[i]][y_low[i] - start_year] for i in range(data_length)])
    z_low_latent = np.array([low_latent[x_low[i]][y_low[i] - start_year] for i in range(data_length)])
    z_low_sign = (np.array(z_low_population) >= 0)
    z_low = np.array([low_percentage[x_low[i]][y_low[i] - start_year] for i in range(data_length)])

    z_difference_infected = list(np.array(z_high_infected) - np.array(z_low_infected))

    ax4 = fig.add_subplot(2, 3, 4, projection='3d')
    ax4.set_xlabel("age")
    ax4.set_ylabel("year")
    ax4.set_zlabel("# latent")

    ax5 = fig.add_subplot(2, 3, 5, projection='3d')
    ax5.set_xlabel("age")
    ax5.set_ylabel("year")
    ax5.set_zlabel("# infected")

    ax6 = fig.add_subplot(2, 3, 6, projection='3d')
    ax6.set_xlabel("age")
    ax6.set_ylabel("year")
    ax6.set_zlabel("infected difference")
    
    ax4.plot_trisurf(x_low, y_low, z_low_latent, cmap=cm.coolwarm, linewidth=0.2, antialiased=True)
    ax5.plot_trisurf(x_low, y_low, z_low_infected, cmap=cm.coolwarm, linewidth=0.2, antialiased=True)
    ax6.plot_trisurf(x_low, y_low, z_difference_infected, cmap=cm.coolwarm, linewidth=0.2, antialiased=True)

    z_difference = z_high - z_low
    z_difference = np.resize(z_difference, (n, yrsOfAnalysis))

    infected_difference = high_infected - low_infected

    population_difference = z_high_population - z_low_population
    population_difference = np.resize(population_difference, (n, yrsOfAnalysis))


    # ========= difference between low and high fidelity ==========
    print("sum of new beta: {0}".format(sum(sum(new_beta))))

    if enable_plot:
        plt.show()

    total_variation = np.sum(np.abs(z_difference[30:60,:]))
    total_infected_difference = np.sum(np.abs(infected_difference))
    total_population_difference = np.sum(np.abs(population_difference[30:60,:]))
    print(file_index)
    print("total variation: {0}".format(total_variation))
    print("infected difference: {0}".format(total_infected_difference))
    print("population difference: {0}".format(total_population_difference))
    print("simulation total infected {0}".format(np.sum(high_infected[:, range(yrsOfAnalysis)])))
    print("approximated total infected {0}".format(np.sum(low_infected)))
    

    total_high_infected = np.sum(high_infected[:, range(yrsOfAnalysis)])
    total_low_infected = np.sum(low_infected[:, range(yrsOfAnalysis)])

    #f_record = open("result/second/test_{0}.csv".format(str(testing_sample)), "a")
    #f_record.write("simulation, {0}, low model, {1}, infected difference (total variation), {2}, population difference, {3}, policy variance, {4}\n".format(total_high_infected, total_low_infected, total_infected_difference, total_population_difference, rerun_variance))


    # ==================== difference analysis ====================
    total_infected_difference_list = []
    total_latent_difference_list = []
    total_population_difference_list = []
    simulation_total_infected_list = []
    approximated_total_infected_list = []
    beta_list = []
    alpha_list = []
    recompute_beta = False

    for sample_i in range(samples_number):
        high_infected = normalized_high_infected_list[sample_i]
        high_population = normalized_high_population_list[sample_i]
        high_latent = normalized_high_latent_list[sample_i]

        I_infected[:, 0] = np.resize(high_infected[:,0], (n))
        N_population[:, 0] = np.resize(high_population[:,0], (n))
        E_latent[:, 0] = np.resize(high_latent[:,0], (n))
        S_safe[:, 0] = N_population[:, 0] - I_infected[:, 0] - E_latent[:,0]

        birth_rate = sfw.zeros((n, yrsOfAnalysis))
        newI = np.zeros((n, yrsOfAnalysis))
        for i in range(yrsOfAnalysis-1):
            birth_rate[0, i] = high_population[0,i+1]
            newI[0,i] = high_infected[0,i+1]

        if recompute_beta:
            previous_samples_number = sample_i
            x_train = np.concatenate([x_train_infected[:,:previous_samples_number * (yrsOfAnalysis - 1),:], x_train_latent[:, :previous_samples_number * (yrsOfAnalysis - 1),:]], axis=1)
            y_train = np.concatenate([y_train_infected[:,:previous_samples_number * (yrsOfAnalysis - 1)], y_train_latent[:, :previous_samples_number * (yrsOfAnalysis - 1)]], axis=1)
    
            weight = np.concatenate( (np.ones((n, previous_samples_number * (yrsOfAnalysis - 1))) ,                 # infected part
                                      np.ones((n, previous_samples_number * (yrsOfAnalysis - 1)))/100  ), axis=1)    # latent part
            beta, alpha = quadratic_solver(x_train, y_train, weight, yrsOfAnalysis)
            beta_list.append(beta)
            alpha_list.append(alpha)



        low_infected, low_population, low_latent = query_low(n, yrsOfAnalysis, G, S_safe, newE, newI, nu, mu, death_rate, alpha, beta, N_population, I_infected, E_latent, birth_rate) # query_low(n, optTimeHorizon, G, S, newE, newI, nu, mu, d, alpha_fast, alpha_slow, beta, N, I, b):

        infected_difference = high_infected - low_infected
        latent_difference = high_latent - low_latent
        population_difference = high_population - low_population

        #total_infected_difference_list.append(np.sum(np.abs(infected_difference)))
        total_infected_difference_list.append(np.sum(np.square(infected_difference)))
        #total_latent_difference_list.append(np.sum(np.abs(latent_difference)))
        total_latent_difference_list.append(np.sum(np.square(latent_difference)))
        #total_population_difference_list.append(np.sum(np.abs(population_difference[30:60,:])))
        total_population_difference_list.append(np.sum(np.square(population_difference[30:60,:])))
        simulation_total_infected_list.append(np.sum(high_infected))
        approximated_total_infected_list.append(np.sum(low_infected))

    print("infected square error: {0}".format(total_infected_difference_list))
    print("latent square error: {0}".format(total_latent_difference_list))
    print("population square error: {0}".format(total_population_difference_list))
    print("simulation total infected {0}".format(simulation_total_infected_list))
    print("approximated total infected {0}".format(approximated_total_infected_list))
    


    #"""
