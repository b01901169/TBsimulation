def diag(a):
    import numpy as np
    return np.matrix(np.diag(a))

def zeros(shape):
    import numpy as np
    return np.matrix(np.zeros(shape))

def ones(shape):
    import numpy as np
    return np.matrix(np.ones(shape))


def load_data(disease, mat_path="SIS_forBryan_TB.mat"):
    '''
    Loads .mat files provided by Sze as numpy matricies/arrays
    '''
    import scipy as sp
    import scipy.io
    import numpy as np
    if disease == 'tb':
        data = sp.io.loadmat(mat_path, struct_as_record=False, squeeze_me=True)
    elif disease == 'gon':
        data = sp.io.loadmat(mat_path, struct_as_record=False, squeeze_me=True)
    else:
        raise Exception('bad disease name: ' + disease)
    beta_array = data['beta']
    beta = []
    for i in range(beta_array.shape[0]):
        beta.append(np.matrix(beta_array[i]))
    f = data['f']
    p = data['p']
    n = data['n']
    G = np.matrix(data['G'])
    for j in range(len(f.I)):
        f.I[j] = np.matrix(f.I[j])
    return f, p, beta, n, G
    
def evaluate(n, T, G, all_S, c, newE, newI, nu, mu, d, alpha_fast, alpha_slow, all_beta, all_N, all_I, b, L):
    import numpy as np
    infected_nu = infected_sim(n, T, G, all_S, c, newE, newI, nu, mu, d, alpha_fast, alpha_slow, all_beta, all_N, all_I, b)
    infected_before = infected_sim(n, T, G, all_S, c, newE, newI, L, mu, d, alpha_fast, alpha_slow, all_beta, all_N, all_I, b)
    return infected_before - infected_nu

def infected(n, T, G, all_S, c, newE, newI, nu, mu, d, alpha_fast, alpha_slow, all_beta, all_N, all_I):
    import numpy as np
#    vals = np.zeros((len(all_beta)))
    for j, (beta, S, N, I) in enumerate(zip(all_beta, all_S, all_N, all_I)):
        transition = [0]*(T+1)
        transition[0] = np.eye(2*n+1)
        x0 = np.transpose(np.bmat([np.transpose(I[:,0]), zeros((n)), ones((1))]))
#        print(x0.sum())
        for t in range(1, T+1):
            block1 = G*diag(1-nu)*diag(1-d) + G*diag(alpha_fast)*diag(S[:,t-1])*diag(1-mu[:,t-1])*beta*diag(np.array(ones((n))/N[:,t-1])[0,:])
            block2 = G*diag(1-mu[:,t-1])*diag(alpha_slow)
            block3 = G*diag(1-alpha_fast)*diag(S[:,t-1])*diag(1-mu[:,t-1])*beta*diag(np.array(ones((n))/N[:,t-1])[0,:])
            block4 = G*diag(1-mu[:,t-1])*diag(1-alpha_slow) 
            allBlocks = np.bmat([[block1, block2, np.transpose(np.matrix(f.newI[:,t-1]))], [block3, block4, np.transpose(np.matrix(newE[:,t-1]))], [zeros((1, n)), zeros((1, n)), ones(1)]])
            transition[t] = allBlocks*transition[t-1];
        return transition
#        print(c.shape, transition[T].shape, x0.shape)
#        vals[j] = (np.transpose(c) * transition[T] * x0)[0,0]
#    return vals


def mean_transition(n, T, G, all_S, c, newE, newI, nu, mu, d, alpha_fast, alpha_slow, all_beta, all_N, all_I):
    import numpy as np
#    vals = np.zeros((len(all_beta)))
    meanTransition = np.zeros((2*n+1, 2*n+1))
    for j, (beta, S, N, I) in enumerate(zip(all_beta, all_S, all_N, all_I)):
        transition = [0]*(T+1)
        transition[0] = np.eye(2*n+1)
        x0 = np.transpose(np.bmat([np.transpose(I[:,0]), zeros((n)), ones((1))]))
#        print(x0.sum())
        for t in range(1, T+1):
            block1 = G*diag(1-nu)*diag(1-d) + G*diag(alpha_fast)*diag(S[:,t-1])*diag(1-mu[:,t-1])*beta*diag(np.array(ones((n))/N[:,t-1])[0,:])
            block2 = G*diag(1-mu[:,t-1])*diag(alpha_slow)
            block3 = G*diag(1-alpha_fast)*diag(S[:,t-1])*diag(1-mu[:,t-1])*beta*diag(np.array(ones((n))/N[:,t-1])[0,:])
            block4 = G*diag(1-mu[:,t-1])*diag(1-alpha_slow) 
            allBlocks = np.bmat([[block1, block2, np.transpose(np.matrix(newI[:,t-1]))], [block3, block4, np.transpose(np.matrix(newE[:,t-1]))], [zeros((1, n)), zeros((1, n)), ones(1)]])
            transition[t] = allBlocks*transition[t-1];
        meanTransition += transition[-1]
    meanTransition /= len(all_beta)
    return meanTransition


def mean_transition_single(n, T, G, all_S, c, newE, newI, nu, mu, d, alpha_fast, alpha_slow, all_beta, all_N, all_I):
    import numpy as np
#    vals = np.zeros((len(all_beta)))
    meanTransition = np.zeros((2*n+1, 2*n+1))
    for j, (beta, S, N, I) in enumerate(zip(all_beta, all_S, all_N, all_I)):
        transition = [0]*(T+1)
        transition[0] = np.eye(2*n+1)
        x0 = np.transpose(np.bmat([np.transpose(I[:,0]), zeros((n)), ones((1))]))
#        print(x0.sum())
        for t in range(1,2):
            block1 = G*diag(1-nu)*diag(1-d) + G*diag(alpha_fast)*diag(S[:,t-1])*diag(1-mu[:,t-1])*beta*diag(np.array(ones((n))/N[:,t-1])[0,:])
            block2 = G*diag(1-mu[:,t-1])*diag(alpha_slow)
            block3 = G*diag(1-alpha_fast)*diag(S[:,t-1])*diag(1-mu[:,t-1])*beta*diag(np.array(ones((n))/N[:,t-1])[0,:])
            block4 = G*diag(1-mu[:,t-1])*diag(1-alpha_slow) 
            allBlocks = np.bmat([[block1, block2, np.transpose(np.matrix(newI[:,t-1]))], [block3, block4, np.transpose(np.matrix(newE[:,t-1]))], [zeros((1, n)), zeros((1, n)), ones(1)]])
            transition[t] = allBlocks*transition[t-1];
        meanTransition += transition[-1]
    meanTransition /= len(all_beta)
    return meanTransition

#        print(c.shape, transition[T].shape, x0.shape)
#        vals[j] = (np.transpose(c) * transition[T] * x0)[0,0]
#    return vals



def spectral(n, T, G, all_S, c, newE, newI, mu, d, alpha_fast, alpha_slow, all_beta, all_N, all_I, L, U, K):
    import numpy as np
    import scipy as sp
    nu = np.copy(L)
    curr = 0
    while (nu - L).sum() < K and curr < n:
        eligible = np.where(nu < U)
        vals = []
        for i in eligible[0]:
            nu_i = nu.copy()
            nu_i[i] += min([U[i] - L[i], K - (nu - L).sum()])
            Bt = mean_transition(n, T, G, all_S, c, newE, newI, nu_i, mu, d, alpha_fast, alpha_slow, all_beta, all_N, all_I)
            vals.append(np.linalg.eigvals(Bt).max())
        to_add = eligible[0][np.argmin(vals)]
        nu[to_add] += min([U[to_add] - L[to_add], K - (nu - L).sum()])
        curr += 1
    return nu

def expect_eig(n, T, G, all_S, c, newE, newI, mu, d, alpha_fast, alpha_slow, all_beta, all_N, all_I, L, U, K):
    import numpy as np
    import scipy as sp
    nu = np.copy(L)
    curr = 0
    while (nu - L).sum() + 0.00001 < K and curr < n:
        ineligible = np.where(nu >= U)
        Bt = mean_transition(n, T, G, all_S, c, newE, newI, nu, mu, d, alpha_fast, alpha_slow, all_beta, all_N, all_I)
        w, v = np.linalg.eig(Bt)
        max_eig = np.argmax(w)
        max_vector = v[:, max_eig]
        max_vector[ineligible] = 0
        to_add = np.argmax(max_vector**2)
        nu[to_add] += min([U[to_add] - L[to_add], K - (nu - L).sum()])
        curr += 1
    return nu
    

def infected_sim(n, T, G, all_S, c, newE, newI, nu, mu, d, alpha_fast, alpha_slow, all_beta, all_N, all_I, b):
    import numpy as np
    vals = np.zeros((len(all_beta), T))
    for j in range(len(all_beta)):
        I_sim = zeros((n, T))
        I_sim[:, 0] = all_I[j][:,0]
        N = zeros((n, T))
        N[:, 0] = np.transpose(np.matrix(all_N[j][:, 0]))
        S = zeros((n, T))
        S[:, 0] = np.transpose(np.matrix(all_S[j][:, 0]))
        beta = all_beta[j]
        vals[j, 0] = I_sim[:, 0].sum()
        for t in range(1, T):
            newExposed = diag(alpha_fast)*diag(np.array(S[:,t-1])[:,0])*diag(1-mu[:,t-1])*beta*diag(np.array(ones((n))/N[:,t-1])[0,:])*I_sim[:,t-1]
            I_sim[:, t] = G*(diag(1-nu)*diag(1-d)*I_sim[:,t-1] + newExposed) + np.transpose(np.matrix(newI[:, t-1]))
            S[:, t] = np.transpose(np.matrix(b[:,t-1])) + G*(diag(1-mu[:,t-1])*S[:,t-1] + diag(nu)*diag(1-d)*I_sim[:,t-1] - newExposed)
            N[:, t] = S[:, t] + I_sim[:, t]
            vals[j, t] = I_sim[:, t].sum()
    return vals



def gradient(n, T, G, S, x0, c, newE, newI, nu, mu, d, alpha_fast, alpha_slow, beta, N):
    import numpy as np
    transition = [0]*(T+1)
    transition[0] = np.eye(2*n+1)
    for t in range(1, T+1):
        block1 = G*diag(1-nu)*diag(1-d) + G*diag(alpha_fast)*diag(S[:,t-1])*diag(1-mu[:,t-1])*beta*diag(np.array(ones((n))/N[:,t-1])[0,:])
        block2 = G*diag(1-mu[:,t-1])*diag(alpha_slow)
        block3 = G*diag(1-alpha_fast)*diag(S[:,t-1])*diag(1-mu[:,t-1])*beta*diag(np.array(ones((n))/N[:,t-1])[0,:])
        block4 = G*diag(1-mu[:,t-1])*diag(1-alpha_slow) 
        allBlocks = np.bmat([[block1, block2, np.transpose(np.matrix(newI[:,t-1]))], [block3, block4, np.transpose(np.matrix(newE[:,t-1]))], [zeros((1, n)), zeros((1, n)), ones(1)]])
#        transition[t] = allBlocks*transition[t-1];
        transition[t] = allBlocks
	
    B_prefix = [0]*(T+1)
    B_prefix[0] = np.eye(2*n+1);
    for i,t in zip(range(1, T+1), range(1, T+1)[::-1]):
    	B_prefix[i] = B_prefix[i-1]*transition[t]
	
    B_suffix= [0]*(T+1)
    B_suffix[0] = np.eye(2*n+1);
    for t in range(1, T+1):
        B_suffix[t] = transition[t]*B_suffix[t-1]
        
    grad = np.zeros((n))
    for i in range(0, n):
        total = zeros((2*n+1, 2*n+1))
        J = zeros((2*n+1, 2*n+1))
        J[i+1, i] = 1
        for t in range(0,T):
            total = total + B_prefix[t] * J * B_suffix[T - t - 1]
        grad[i] = (1 - d[i])*np.trace((c * np.transpose(x0)) * total);
    return grad 

def gradient_sum(n, T, G, S, x0, c, newE, newI, nu, mu, d, alpha_fast, alpha_slow, beta, N):
    import numpy as np
    transition = [0]*(T+1)
    transition[0] = np.eye(2*n+1)
    for t in range(1, T+1):
        block1 = G*diag(1-nu)*diag(1-d) + G*diag(alpha_fast)*diag(S[:,t-1])*diag(1-mu[:,t-1])*beta*diag(np.array(ones((n))/N[:,t-1])[0,:])
        block2 = G*diag(1-mu[:,t-1])*diag(alpha_slow)
        block3 = G*diag(1-alpha_fast)*diag(S[:,t-1])*diag(1-mu[:,t-1])*beta*diag(np.array(ones((n))/N[:,t-1])[0,:])
        block4 = G*diag(1-mu[:,t-1])*diag(1-alpha_slow) 
        allBlocks = np.bmat([[block1, block2, np.transpose(np.matrix(newI[:,t-1]))], [block3, block4, np.transpose(np.matrix(newE[:,t-1]))], [zeros((1, n)), zeros((1, n)), ones(1)]])
#        transition[t] = allBlocks*transition[t-1];
        transition[t] = allBlocks
	
    B_prefix = [0]*(T+1)
    B_prefix[0] = np.eye(2*n+1);
    for i,t in zip(range(1, T+1), range(1, T+1)[::-1]):
    	B_prefix[i] = B_prefix[i-1]*transition[t]
	
    B_suffix= [0]*(T+1)
    B_suffix[0] = np.eye(2*n+1);
    for t in range(1, T+1):
        B_suffix[t] = transition[t]*B_suffix[t-1]
        
    grad = np.zeros((n))
    for i in range(0, n):
        total = zeros((2*n+1, 2*n+1))
        J = zeros((2*n+1, 2*n+1))
        J[i+1, i] = 1
        for Tprime in range(1, T+1):
            for t in range(0,Tprime):
                total = total + B_prefix[(T - Tprime) + t] * J * B_suffix[T - t - 1]
            grad[i] += (1 - d[i])*np.trace((c * np.transpose(x0)) * total);
    return grad 


def greedy(grad, U, L, K):
    '''
    Greedily select budget number of elements with highest weight according to
    grad
    '''
    import numpy as np
    elements = range(len(grad))
    combined = zip(elements, grad)
    combined.sort(key = lambda x: -x[1])
    sorted_groups = [x[0] for x in combined]
    nu = np.copy(L)
    curr = 0
    while (nu - L).sum() < K and curr < len(grad):
        amount_add = min([U[sorted_groups[curr]] - L[sorted_groups[curr]], K - (nu - L).sum()])
        nu[[sorted_groups[curr]]] += amount_add
        curr += 1
    return nu


def stochastic_frank_wolfe(n, T, G, S, c, newE, newI, mu, d, alpha_fast, alpha_slow, beta, N, L, U, K, num_iter, I):
    import numpy as np
    import random
    nu = np.zeros((n))
    print "total iterations: {0}".format((num_iter, len(beta)))
    for i in range(num_iter):
        grad = np.zeros((n))
#        for j in range(len(beta)):
        for j in random.sample(range(len(beta)), 100):
            print (i, j)
            x0 = np.transpose(np.bmat([np.transpose(I[j][:,0]), zeros((n)), ones((1))]))
#            grad = gradient(n, T, G, S[j], x0, c, newE, newI, nu+L, mu, d, alpha_fast, alpha_slow, beta[j], N[j])
            #Swap this line in for Tth step only, vs sum of 1...T 
#            grad = grad + gradient(n, T, G, S[j], x0, c, newE, newI, nu+L, mu, d, alpha_fast, alpha_slow, beta[j], N[j])
            grad = grad + gradient_sum(n, T, G, S[j], x0, c, newE, newI, nu+L, mu, d, alpha_fast, alpha_slow, beta[j], N[j])
        grad = grad / len(beta)
        v = greedy(grad, U - L, np.zeros((n)), K)
        nu = nu + (1./num_iter)*v
    nu = nu + L;
    return nu
