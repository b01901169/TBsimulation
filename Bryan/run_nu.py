import numpy as np
import sys
from sfw import load_data, stochastic_frank_wolfe, greedy, ones, zeros, evaluate
import itertools
import multiprocessing
from functools import partial
import pickle

def run_param(param, f, p, beta, n, G, L):
    optTimeHorizon, K = param
    num_iter = 100
    c = np.bmat([[ones((n, 1))], [zeros((n+1, 1))]])
    U = np.ones((n)) * 0.05;  
    U = L + U
    #greedy on highest prevalence
    nus = {}
    nuPrev = greedy(np.array(f.I[0][:, 0])[:,0], U, L, K)
    nus['prev'] = nuPrev
    #equal split of budget across groups
    nuEqual = L + np.ones((n)) * K/n
    nus['equal'] = nuEqual
    #greedy on highest degree in beta matrix
    meanBeta = beta[0]
    for x in beta[1:]:
        meanBeta += x
    meanBeta /= len(beta)
    degrees = meanBeta.sum(axis=1)
    nuDegree = greedy(degrees, U, L, K)
    nus['degree'] = nuDegree
    #SFW-based
    print "running SFW based method ..."
    nuSFW = stochastic_frank_wolfe(n, optTimeHorizon, G, f.S, c, f.newE, f.newI, f.mu, f.d, f.alpha_fast, f.alpha_slow, beta, f.N, L, U, K, num_iter, f.I)
    nuMyopic = stochastic_frank_wolfe(n, 1, G, f.S, c, f.newE, f.newI, f.mu, f.d, f.alpha_fast, f.alpha_slow, beta, f.N, L, U, K, 1, f.I)
    nus['sfw'] = nuSFW
    nus['myopic'] = nuMyopic
    print "running evaluation ..."
    #evaluate all of the nus
    algos = ['prev', 'equal', 'degree', 'sfw', 'myopic']
    averted = {}
    for name in algos:
        averted[name] = evaluate(n, optTimeHorizon, G, f.S, c, f.newE, f.newI, nus[name], f.mu, f.d, f.alpha_fast, f.alpha_slow, beta, f.N, f.I, f.b, L)
    return nus, averted


if __name__ == '__main__':
#    disease = sys.argv[1]
    disease = 'tb'
    print "loading ..."
    f, p, beta, n, G = load_data(disease)
    print "f: {0}".format(f)
    print "p: {0}".format(p)
    print "beta: {0}".format(len(beta))
    print "n: {0}".format(n)
    print "G: {0}".format(G)
    print "finished loading ..."
   
    if disease == 'gon':
        L = p.nu_sampled[1][:, p.T-1]
    elif disease == 'tb':
        L = p.nu[:,p.T-1]
    else:
        raise Exception('bad disease name')
    Ks = np.linspace(0.05, 0.4, 10)
    #Ks = [Ks[1], Ks[7]]
    Ks = [Ks[7]]
    Ts = [25]
#    params = itertools.product(range(1, f.T+1), np.linspace(0.05, 0.4, 10))
    params = itertools.product(Ts, Ks)

    print "p.T: {0}".format(p.T)
    print "L: {0}".format(L)

'''
    print "start computing ..."
    pool = multiprocessing.Pool(1)
    rp = partial(run_param, f=f, p=p, beta=beta, n=n, G=G, L=L)
    results = pool.map(rp, params)
    print "finished computing ..."
    
#    with open('/home/rcf-proj2/mj1/bwilder/results_sum_small_' + disease, 'wb') as f:
    with open('results_sum_small_' + disease, 'wb') as f:
        pickle.dump(results, f)

#'''
