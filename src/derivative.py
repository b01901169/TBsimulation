import numpy as np


def getDerivativeByNu(n, yrsOfAnalysis, I, S, newI, birth_rate, death_rate, nu, mu, beta):
    derI = np.zeros((n, n, yrsOfAnalysis))
    derS = np.zeros((n, n, yrsOfAnalysis))

    #for j in range(n): # computing the derivative of nu[j]
    for t in range(yrsOfAnalysis-1):
        for i in range(n-1):
            #print t,i
            newExposed = np.array([beta[i][k] * (derS[:,k,t] * I[k][t]**2 + derI[:,k,t] * S[k][t]**2) / float((I[k][t] + S[k][t])**2)  for k in range(n)])
            newExposed[np.isnan(newExposed)] = 0
            newExposed = np.sum(newExposed, axis=0)

            derI[:,i+1,t+1] = derI[:,i,t] * (1 - death_rate[i]) * (1 - nu[i]) + (1 - mu[i][t]) * newExposed
            derS[:,i+1,t+1] = derS[:,i,t] * (1 - mu[i][t]) + derI[:,i,t] * (1 - death_rate[i]) * nu[i] - (1 - mu[i][t]) * newExposed

            derI[i][i+1][t+1] += - I[i][t] * (1 - death_rate[i])
            derS[i][i+1][t+1] += I[i][t] * (1 - death_rate[i])

    return derI, derS

def getDerivativeByInfected(n, yrsOfAnalysis, I, S, newI, birth_rate, death_rate, nu, mu, beta):
    derI = np.zeros((n, n, yrsOfAnalysis))
    derS = np.zeros((n, n, yrsOfAnalysis))
    for i in range(n):
        derI[i,i,0] = float(1)

    #for j in range(n): # computing the derivative of I[:,0]
    for t in range(yrsOfAnalysis-1):
        for i in range(n-1):
            #print t,i
            newExposed = np.array([beta[i][k] * (derS[:,k,t] * I[k][t]**2 + derI[:,k,t] * S[k][t]**2) / float((I[k][t] + S[k][t])**2)  for k in range(n)])
            newExposed[np.isnan(newExposed)] = 0
            newExposed = np.sum(newExposed, axis=0)

            derI[:,i+1,t+1] = derI[:,i,t] * (1 - death_rate[i]) * (1 - nu[i]) + (1 - mu[i][t]) * newExposed
            derS[:,i+1,t+1] = derS[:,i,t] * (1 - mu[i][t]) + derI[:,i,t] * (1 - death_rate[i]) * nu[i] - (1 - mu[i][t]) * newExposed

    return derI, derS

def getDerivativeBySafe(n, yrsOfAnalysis, I, S, newI, birth_rate, death_rate, nu, mu, beta):
    derI = np.zeros((n, n, yrsOfAnalysis))
    derS = np.zeros((n, n, yrsOfAnalysis))
    for i in range(n):
        derS[i,i,0] = 1

    #for j in range(n): # computing the derivative of S[:,0]
    for t in range(yrsOfAnalysis-1):
        for i in range(n-1):
            #print t,i
            newExposed = np.array([beta[i][k] * (derS[:,k,t] * I[k][t]**2 + derI[:,k,t] * S[k][t]**2) / float((I[k][t] + S[k][t])**2)  for k in range(n)])
            newExposed[np.isnan(newExposed)] = 0
            newExposed = np.sum(newExposed, axis=0)

            derI[:,i+1,t+1] = derI[:,i,t] * (1 - death_rate[i]) * (1 - nu[i]) + (1 - mu[i][t]) * newExposed
            derS[:,i+1,t+1] = derS[:,i,t] * (1 - mu[i][t]) + derI[:,i,t] * (1 - death_rate[i]) * nu[i] - (1 - mu[i][t]) * newExposed

    return derI, derS
