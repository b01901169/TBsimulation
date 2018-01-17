import numpy as np
import itertools
import sys

Bryan_path = "/home/kai/Dropbox/USC/publication/TBsimulationCode/Bryan/"
sys.path.insert(0, Bryan_path)

import sfw

if __name__ == "__main__":
    disease = 'tb'
    print "loading ..."
    f, p, beta, n, G = sfw.load_data(disease, mat_path=Bryan_path + "SIS_forBryan_TB.mat")
    c = np.bmat([[sfw.ones((n, 1))], [sfw.zeros((n+1, 1))]])

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
    K = Ks[7]
    optTimeHorizon = 25
#    params = itertools.product(range(1, f.T+1), np.linspace(0.05, 0.4, 10))

    print "p.T: {0}".format(p.T)
    print "L: {0}".format(L)
    nu = L

    infected_nu = sfw.infected_sim(n, optTimeHorizon, G, f.S, c, f.newE, f.newI, nu, f.mu, f.d, f.alpha_fast, f.alpha_slow, beta, f.N, f.I, f.b)

    print infected_nu
