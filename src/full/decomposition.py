import numpy as np
import scipy

from matplotlib import pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import (RBF, Matern, RationalQuadratic,
                                              ExpSineSquared, DotProduct,
                                              ConstantKernel, WhiteKernel)



class Decomposition:
    def __init__(self, J, fList, g, kernelList=None):
        # J: number of decomposed subfunctions
        # fList: list (size J) of subfunctions. The function is either a lambda function or a callable function with access to f(x)
        self.J = J
        self.fList = fList
        self.g = g
        self.kernelList = kernelList

        
    def get_function_value(self, x):
        # input: x (d dimension array)
        # output: function value f(x) = g(f1(x), f2(x), ..., fJ(x))
        subfunction_values = [0] * self.J
        for i in range(self.J):
            subfunction_values[i] = fList[j](x)

        g_value = self.g(subfunction_values)
        return g_value


    def get_subfunction_value(self, x):
        # input: x (d dimension array)
        # output: J values representing the values of individual subfunctions
        subfunction_values = [0] * self.J
        for i in range(self.J):
            subfunction_values[i] = fList[j](x)

        return subfunction_values


class GPUCB:
    def __init__(self, f, kernel, dimension, upper_bound, constraints, delta=0.05, a=1, b=1, gp_alpha=0.01, initial_point=0.3):
        self.f = f
        self.kernel = kernel
        self.T = 1
        self.dimension = dimension
        self.upper_bound = upper_bound
        self.constraints = constraints
        self.delta = delta
        self.a = a
        self.b = b
        self.gp_alpha = gp_alpha

        self.sample_points = np.reshape([initial_point], (1, dimension))
        self.sample_values = np.reshape([f(initial_point)], (1,1))
        self.bds = [(0, upper_bound) for i in range(dimension)]

    def run(self, iterations):
        for iteration in range(iterations):
            beta_t = 2 * np.log(2 * self.T**2 * np.pi**2 / (3 * self.delta)) + 2 * self.dimension * np.log(self.T**2 * self.dimension * self.b * self.upper_bound * np.sqrt(np.log(4 * self.dimension * self.a / self.delta)))

            gpr = GaussianProcessRegressor(kernel=self.kernel, optimizer=None, alpha=self.gp_alpha)
            gpr.fit(self.sample_points, self.sample_values)

            fn = lambda x: -self.GPUCB_objective_value(gpr, beta_t, x)
            res = scipy.optimize.minimize(fn, self.sample_points[-1], bounds=self.bds, constraints=self.constraints)

            print res
            new_x = res.x
            new_fun = res.fun
            new_objective_value = f(new_x)
            print ("iteration: {0}, new sample point: {1}, objective value: {2}, beta: {3}".format(iteration, new_x, new_objective_value, beta_t))

            self.sample_points = np.concatenate((self.sample_points, np.reshape(new_x, (1, self.dimension))))
            self.sample_values = np.concatenate((self.sample_values, np.reshape(new_objective_value, (1,1))))
            self.T = self.T + 1

    def predict(self, x):
        gpr = GaussianProcessRegressor(kernel=self.kernel, optimizer=None, alpha=self.gp_alpha)
        gpr.fit(self.sample_points, self.sample_values)

        mean, std = gpr.predict(np.reshape(x, (len(x), self.dimension)), return_std=True)
        return mean, std

    def get_beta_t(self):
        beta_t = 2 * np.log(2 * self.T**2 * np.pi**2 / (3 * self.delta)) + 2 * self.dimension * np.log(self.T**2 * self.dimension * self.b * self.upper_bound * np.sqrt(np.log(4 * self.dimension * self.a / self.delta)))
        return beta_t


    def GPUCB_objective_value(self, gpr, beta_t, x):
        mean, std = gpr.predict(np.reshape(x, (len(x), self.dimension)), return_std=True)
        return mean + np.sqrt(beta_t) * std



if __name__ == "__main__":
    upper_bound = 1
    grid_size = 1000
    kernel = RBF(length_scale=0.3)
    gp_alpha = 0.01
    dimension = 1
    random_seed = np.random.randint(1000)
    #budget_constraint = 0.3
    #constraints = ({"type": "eq", "fun": lambda x: np.sum(x) - budget_constraint})
    constraints = None

    X_ = np.reshape(np.linspace(0, upper_bound, grid_size), (grid_size, 1))

    gp = GaussianProcessRegressor(kernel=kernel, optimizer=None, alpha=gp_alpha)
    target = gp.sample_y(X_, 1, random_seed)

    def f(x):
        index = np.argmax(x < X_)
        value = (index / float(grid_size) - x) * target[index-1] + (1 - index / float(grid_size) + x) * target[index]
        return value

    GPUCBsolver = GPUCB(f, kernel, dimension, upper_bound, constraints)
    #GPUCBsolver.run(1)




