from sklearn.gaussian_process.kernels import Kernel, _approx_fprime, Hyperparameter
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import pairwise_kernels

from sklearn.gaussian_process.kernels import RBF, Matern

class composedKernel(Kernel):
    def __init__(self, coefficients, subkernels): # respectively n, T, g_i(x), k_i(x,x')
        assert(len(coefficients) == len(subkernels))

        self.time_horizon = len(coefficients)
        self.coefficients = coefficients
        self.subkernels = subkernels

    def __call__(self, X, Y=None): # eval_gradient=False # not implemented yet
        X = np.atleast_2d(X)
        if Y is None:
            K = np.zeros((len(X), len(X)))
        else:
            K = np.zeros((len(X), len(Y)))
        for i in range(self.time_horizon):
            if Y is None:
                subK = np.matrix(self.subkernels[i](X))
                DX = np.matrix(np.diag(self.coefficients[i](X)))
                DY = np.matrix(np.diag(self.coefficients[i](X)))
            else:
                subK = np.matrix(self.subkernels[i](X, Y))
                DX = np.matrix(np.diag(self.coefficients[i](X)))
                DY = np.matrix(np.diag(self.coefficients[i](Y)))
            K += np.array(DX * subK * DY)
        return K

    def is_stationary(self):
        return False

    def diag(self, X):
        total_diag = np.zeros(len(X))
        for i in range(self.time_horizon):
            sub_diag = self.subkernels[i].diag(X)
            coeff = self.coefficients[i](X)
            total_diag += coeff * sub_diag * coeff
        return total_diag

if __name__ == "__main__":
    f1 = lambda X: [2 for i in range(len(X))]
    f2 = lambda X: [1 for i in range(len(X))]
    f3 = lambda X: [3 for i in range(len(X))]
    k1 = RBF()
    k2 = RBF(length_scale=2.0)
    k3 = Matern()

    coefficients = [f1, f2, f3]
    subkernels = [k1, k2, k3]

    X = np.array([[1,2,3,4]])
    Y = np.array([[4,3,2,1]])

    ck = composedKernel(3, coefficients, subkernels)
