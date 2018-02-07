import matlab.engine
import os
import uptakePolicy
import numpy as np
import itertools
import pickle
import sys
import argparse

# =============== matplot ================
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

# ============= our library ==============
Bryan_path = "/home/kai/Dropbox/USC/publication/TBsimulationCode/Bryan/"
dir_path = "/home/kai/Dropbox/USC/publication/TBsimulationCode/"
sys.path.insert(0, Bryan_path)

import sfw

# ============== functions ===============
def query_low(n, optTimeHorizon, G, S, newE, newI, nu, mu, d, alpha_fast, alpha_slow, beta, N, I, b):
    disease = 'tb'
    print "loading ..."
    #f, p, beta, n, G = sfw.load_data(disease, mat_path=Bryan_path + "SIS_forBryan_TB.mat")
    c = np.bmat([[sfw.ones((n, 1))], [sfw.zeros((n+1, 1))]])

    print "f: {0}".format(f)
    print "p: {0}".format(p)
    print "beta: {0}".format(len(beta))
    print "n: {0}".format(n)
    print "G: {0}".format(G)
    print "finished loading ..."

    #if disease == 'gon':
    #    L = p.nu_sampled[1][:, p.T-1]
    #elif disease == 'tb':
    #    L = p.nu[:,p.T-1]
    #else:
    #    raise Exception('bad disease name')
    Ks = np.linspace(0.05, 0.4, 10)
    #Ks = [Ks[1], Ks[7]]
    K = Ks[7]
    optTimeHorizon = 25
#    params = itertools.product(range(1, f.T+1), np.linspace(0.05, 0.4, 10))

    #print "p.T: {0}".format(p.T)
    #print "f.N: {0}".format(f.N)
    #print "L: {0}".format(L)
    #nu = L

    #infected_nu = sfw.infected_sim(n, optTimeHorizon, G, f.S, c, f.newE, f.newI, nu, f.mu, f.d, f.alpha_fast, f.alpha_slow, beta, f.N, f.I, f.b)
    #infected_nu, N_sim = sfw.single_infected_sim(n, optTimeHorizon, G, f.S[0], c, f.newE, f.newI, nu, f.mu, f.d, f.alpha_fast, f.alpha_slow, beta[0], f.N[0], f.I[0], f.b)
    infected_nu, N_sim = sfw.single_infected_sim(n, optTimeHorizon, G, S, c, newE, newI, nu, mu, d, alpha_fast, alpha_slow, beta, N, I, b)

    #print infected_nu, N_sim
    return infected_nu, N_sim



def query_high(NumPpl, yrsOfAnalysis, uptakeUrbanKnowledge, uptakeAgeBracs):
    eng_path = "/home/kai/Dropbox/USC/publication/TBsimulationCode/TBsimulationCode"
    eng = matlab.engine.start_matlab()
    eng.cd(eng_path)
    
    Group1, age_population = eng.TBsimulation_jan23('.','p01', 130+yrsOfAnalysis, NumPpl, '-r0','loadBurnIn', 2018,0, 1.140, 0.00018, 0.00007, 0.0025, 0.8, 1.3,'base', 'customized',
                                    matlab.double(uptakeUrbanKnowledge), matlab.int8(uptakeAgeBracs), nargout=2)

    return Group1, age_population


if __name__ == "__main__":
    disease = 'tb'
    yrsOfAnalysis = 25

    parser = argparse.ArgumentParser(description="filename parser")
    parser.add_argument("-d", "--date", help="Input date")
    parser.add_argument("-m", "--method", help="Input one of the methods: initial, optimal, full, peak40, peak55")

    args = parser.parse_args()
    print "args: {0}".format(args)

    date = args.date
    file_index = args.method

    print "date: {0}, method {1}".format(date, file_index)

    print "loading ..."
    f, p, beta, n, G = sfw.load_data(disease, mat_path=Bryan_path + "SIS_forBryan_TB.mat")

    # =================== parameter setting ===================
    n = 86
    beta = np.matrix(np.ones((n,n))/120)
    #beta = np.matrix(np.zeros((n,n)))
    #beta = beta[0]
    newE = np.zeros((n, yrsOfAnalysis))
    newI = np.zeros((n, yrsOfAnalysis))
    alpha_fast = np.ones(n)
    alpha_slow = np.ones(n)

    G = np.zeros((n,n))
    for i in range(n-1):
        G[i+1][i] = 1
    G = np.matrix(G)
    death_rate = np.array([f.d[0]] * n)
    
    mu = sfw.zeros((n, yrsOfAnalysis))
    I_infected = sfw.zeros((n, yrsOfAnalysis))
    N_population = sfw.zeros((n, yrsOfAnalysis))
    S_safe = sfw.zeros((n, yrsOfAnalysis))
    birth_rate = sfw.zeros((n, yrsOfAnalysis))
    initial_nu = np.zeros(n)
    for i in range(n):
        modified_index = min(max(i-30,0), 30-1)
        modified_birth_index = min(max(i-30,2), 30-1)
        #modified_index = i
        mu[i] = f.mu[modified_index]
        I_infected[i] = f.I[0][modified_index]
        N_population[i] = f.N[0][modified_index]
        S_safe[i] = f.S[0][modified_index]
        birth_rate[i] = f.b[modified_birth_index]
        initial_nu[i] = p.nu[modified_index,p.T-1]
    birth_rate[0] = f.b[0]
        
    d = np.array([f.d[0]] * n)
    #d = np.ones((n))/20

    # ===================== plot setting ======================

    fig = plt.figure()
    # =================== given parameters ====================
    # ------------------ high fidelity model ------------------
    NumPpl = 200000
    yrsOfAnalysis = 25

    # ------------------- low fidelity model ------------------
    N = f.N

    # ================== varied parameters ====================
    I = f.I
    #initial_nu = np.array(p.nu[:,p.T-1]) # TODO
    optimal_nu = np.array(
        [ 0.13294343,  0.13041813,  0.15050087,  0.12541041,  0.12393383, # 5
        0.13294343,  0.13041813,  0.15050087,  0.12541041,  0.12393383, # 10
        0.13294343,  0.13041813,  0.15050087,  0.12541041,  0.12393383, # 15
        0.13294343,  0.13041813,  0.15050087,  0.12541041,  0.12393383, # 20
        0.13294343,  0.13041813,  0.15050087,  0.12541041,  0.12393383, # 25
        0.13294343,  0.13041813,  0.15050087,  0.12541041,  0.12393383, # 30
        0.13294343,  0.13041813,  0.15050087,  0.12541041,  0.12393383, # 35
        0.1225361 ,  0.12121306,  0.11996091,  0.12130266,  0.12221708, # 40
        0.08431062,  0.08534621,  0.08660821,  0.08617747,  0.08530207, # 45
        0.08464646,  0.08420093,  0.08395899,  0.08456231,  0.08522047, # 50
        0.07465641,  0.07582088,  0.07732078,  0.13136248,  0.13945654, # 55
        0.14882135,  0.15973757,  0.17256779,  0.17394538,  0.11919841, # 60
        0.14882135,  0.15973757,  0.17256779,  0.17394538,  0.11919841, # 65
        0.14882135,  0.15973757,  0.17256779,  0.17394538,  0.11919841, # 70
        0.14882135,  0.15973757,  0.17256779,  0.17394538,  0.11919841, # 75
        0.14882135,  0.15973757,  0.17256779,  0.17394538,  0.11919841, # 80
        0.14882135,  0.15973757,  0.17256779,  0.17394538,  0.11919841, 0.11919841] # 86
    )
    #optimal_nu = optimal_nu[30:60]
    peak55_nu = np.zeros(n)
    peak55_nu[26:29] += 1

    peak45_nu = np.zeros(n)
    peak45_nu[15:18] += 1

    peak40_nu = np.zeros(n)
    peak40_nu[11:14] += 1

    peak35_nu = np.zeros(n)
    peak35_nu[4:7] += 1

    full_nu = np.ones(n)
    null_nu = np.zeros(n)

    if file_index == "initial":
        nu = initial_nu
    elif file_index == "optimal":
        nu = optimal_nu
    elif file_index == "peak35":
        nu = peak35_nu
    elif file_index == "peak40":
        nu = peak40_nu
    elif file_index == "peak45":
        nu = peak45_nu
    elif file_index == "peak55":
        nu = peak55_nu
    elif file_index == "full":
        nu = full_nu
    elif file_index == "null":
        nu = null_nu
    else:
        print "method not found"
        raise


    difference_nu = nu - initial_nu
    print "difference sum: {0}, difference: {1}".format(np.sum(difference_nu), difference_nu)

    #"""
    # =================== high fidelity model ==================
    # print "============================ high fidelity ============================="
    # uptakeUrbanKnowledge = np.array(uptakePolicy.finerUptakeUrbanKnowledge)
    # #uptakeUrbanKnowledge[30:60, 0] = nu
    # uptakeUrbanKnowledge[30:60, 0] = nu.copy()
    # #uptakeUrbanKnowledge[30:60, 1] = nu
    # uptakeUrbanKnowledge[30:60, 1] = nu.copy()
    # high_infected, high_population = query_high(NumPpl, yrsOfAnalysis, uptakeUrbanKnowledge.tolist(), uptakePolicy.finerUptakeAgeBracs)
    # pickle.dump((high_infected, high_population), open(dir_path + "data/high_fidelity_{0}_{1}.data".format(date, file_index), "wb"))

    # # ------------------- ploting 3D surface ------------------
    # high_infected = np.array(high_infected)
    # high_population = np.array(high_population)
    # high_percentage = high_infected / high_population

    # x_length, y_length = high_infected.shape
    # x_start_year = 30 # 30 years old started
    # x_end_year = x_length - 50 # 60 years old
    # y_year_length = y_length / 12
    # start_year = 1996 - 130
    # x1 = [i for i in range(x_start_year, x_end_year) for j in range(130, y_year_length)]
    # y1 = [start_year + j for i in range(x_start_year, x_end_year) for j in range(130, y_year_length)]
    # data_length = (x_end_year - x_start_year) * (y_year_length - 130)
    # z1_infected = [high_infected[x1[i]][(y1[i] - start_year) * 12] for i in range(data_length)]
    # z1_population = [high_population[x1[i]][(y1[i] - start_year) * 12] for i in range(data_length)]
    # z1 = [high_percentage[x1[i]][(y1[i] - start_year) * 12] for i in range(data_length)]

    # ax1 = fig.add_subplot(1, 3, 1, projection='3d')
    # ax1.set_xlabel("age")
    # ax1.set_ylabel("year")
    # ax1.set_zlabel("# infected")
    # 
    # ax1.plot_trisurf(x1, y1, z1, cmap=cm.coolwarm, linewidth=0.2, antialiased=True)
    # #plt.show()


    #"""
    #"""
    # =================== low fidelity model ==================
    print "============================ low fidelity ============================="
    print beta
    low_infected, low_population = query_low(n, yrsOfAnalysis, G, S_safe, newE, newI, nu, mu, death_rate, alpha_fast, alpha_slow, beta, N_population, I_infected, birth_rate) # query_low(n, optTimeHorizon, G, S, newE, newI, nu, mu, d, alpha_fast, alpha_slow, beta, N, I, b):
    #infected_nu, N_sim = sfw.single_infected_sim(n, optTimeHorizon, G, f.S[0], c, f.newE, f.newI, nu, f.mu, f.d, f.alpha_fast, f.alpha_slow, beta[0], f.N[0], f.I[0], f.b)
    pickle.dump((low_infected, low_population), open(dir_path + "data/low_fidelity_{0}_{1}.data".format(date, file_index), "wb"))

    # ------------------- ploting 3D surface ------------------
    low_infected = np.array(low_infected)
    low_population = np.array(low_population)
    low_percentage = low_infected / low_population

    x_length, y_length = low_infected.shape
    start_age = 0
    start_year = 1995
    x2 = [i + start_age for i in range(x_length) for j in range(y_length)]
    y2 = [j + start_year for i in range(x_length) for j in range(y_length)]
    data_length = x_length * y_length
    low_infected = np.asarray(low_infected)
    z2_infected = [low_infected[x2[i] - start_age][y2[i] - start_year] for i in range(data_length)]
    z2_population = [low_population[x2[i] - start_age][y2[i] - start_year] for i in range(data_length)]
    z2_sign = (np.array(z2_population) >= 0)
    z2 = [low_percentage[x2[i] - start_age][y2[i] - start_year] for i in range(data_length)]

    ax2 = fig.add_subplot(1, 3, 2, projection='3d')
    ax2.set_xlabel("age")
    ax2.set_ylabel("year")
    ax2.set_zlabel("# infected")
    
    #ax2.plot_trisurf(x2, y2, z2, cmap=cm.coolwarm, linewidth=0.2, antialiased=True)
    #ax2.plot_trisurf(x2, y2, z2_infected, cmap=cm.coolwarm, linewidth=0.2, antialiased=True)
    ax2.plot_trisurf(x2, y2, z2_population, cmap=cm.coolwarm, linewidth=0.2, antialiased=True)
    #ax2.plot_trisurf(x2, y2, z2_sign, cmap=cm.coolwarm, linewidth=0.2, antialiased=True)

    # ========= difference between low and high fidelity ==========
    # z3 = np.array(z2) - np.array(z1)

    # ax3 = fig.add_subplot(1, 3, 3, projection='3d')
    # ax3.set_xlabel("age")
    # ax3.set_ylabel("year")
    # ax3.set_zlabel("# infected")
    # ax3.plot_trisurf(x1, y1, z3, cmap=cm.coolwarm, linewidth=0.2, antialiased=True)

    plt.show()

    # total_variation = np.sum(np.abs(z3))
    print file_index
    print "total variation: {0}".format(total_variation)
    print "simulation total infected {0}".format(np.sum(high_infected[30:60, range(130*12, 155*12, 12)]))
    print "approximated total infected {0}".format(np.sum(low_infected, 1))
    print "approximated total infected {0}".format(np.sum(low_infected))
    

    #"""
