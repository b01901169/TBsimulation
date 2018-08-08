import numpy as np
import scipy

from matplotlib import pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import (RBF, Matern, RationalQuadratic,
                                              ExpSineSquared, DotProduct,
                                              ConstantKernel, WhiteKernel)



class Decomposition:
    def __init__(self, J, fList, g, gList, kernelList):
        # J: number of decomposed subfunctions
        # fList: list (size J) of subfunctions. The function is either a lambda function or a callable function with access to f(x)
        self.J = J
        self.fList = fList
        self.g = g
        self.gList = gList
        self.kernelList = kernelList

        
    def get_function_value(self, x):
        # input: x (d dimension array)
        # output: function value f(x) = g(f1(x), f2(x), ..., fJ(x))
        assert(type(x) == float or type(x) == int or len(x) == 1)

        subfunction_values = np.zeros(self.J)
        for i in range(self.J):
            subfunction_values[i] = self.fList[i](x)

        g_value = self.g(subfunction_values)
        return g_value


    def get_subfunction_values(self, x):
        # input: x (d dimension array)
        # output: J values representing the values of individual subfunctions
        assert(type(x) == float or type(x) == int or len(x) == 1)

        subfunction_values = np.zeros(self.J)
        for i in range(self.J):
            subfunction_values[i] = self.fList[i](x)

        return subfunction_values

    def get_coefficients(self, x):
        assert(type(x) == float or type(x) == int or len(x) == 1)
        coefficients = np.zeros((self.J, len(x)))
        for i in range(self.J):
            coefficients[i] = self.gList[i](x)

        return coefficients

class GPUCB:
    def __init__(self, f, kernel, dimension, upper_bound, constraints, delta=0.05, a=1, b=1, gp_alpha=0.01, initial_point=None, X_=None, discrete=False, linear=True, B=10):
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
        self.discrete = discrete
        self.linear = linear
        self.B = B
        if np.shape(X_):
            grid_size = 1000
            self.X_ = np.reshape(np.linspace(0, upper_bound, grid_size), (grid_size, 1))
        else:
            self.X_ = X_
        self.gpr = None
        self.regret = 0

        # ---------------------- initial sample -----------------------
        if not initial_point and dimension == 1:
            initial_point = np.random.rand() * upper_bound
        else:
            raise ("No specified initial point!")

        self.sample_points = np.reshape([initial_point], (1, dimension))
        self.sample_values = np.reshape([f(initial_point)], (1,1))
        self.bds = [(0, upper_bound) for i in range(dimension)]
        if dimension == 1:
            self.true_optimal = np.max([f(x) for x in self.X_])
            print ("True optimal: {0}".format(self.true_optimal))
        else:
            self.true_optimal = None

    def run(self, iterations):
        for iteration in range(iterations):
            # beta_t = 2 * np.log(2 * self.T**2 * np.pi**2 / (3 * self.delta)) + 2 * self.dimension * np.log(self.T**2 * self.dimension * self.b * self.upper_bound * np.sqrt(np.log(4 * self.dimension * self.a / self.delta)))
            beta_t = self.get_beta_t()

            gpr = GaussianProcessRegressor(kernel=self.kernel, optimizer=None, alpha=self.gp_alpha)
            self.gpr = gpr
            gpr.fit(self.sample_points, self.sample_values)

            fn = lambda x: -self.GPUCB_objective_value(gpr, beta_t, x)
            # ------------------- optimization choice ---------------------
            if self.discrete:
                objective_values = fn(self.X_)
                optimal_index = np.argmin(objective_values)
                optimal_x = self.X_[optimal_index]
                optimal_fun = objective_values[optimal_index]
                optimal_objective_value = self.f(optimal_x)

                new_x = optimal_x
                new_fun = optimal_fun
            else:
                # initial_point = self.sample_points[np.argmax(self.sample_values)]
                initial_point = np.random.rand() * self.upper_bound

                res = scipy.optimize.minimize(fn, initial_point, bounds=self.bds, constraints=self.constraints)
                # print (res)

                new_x = res.x
                new_fun = res.fun

            # ------------------- suboptimal choice -----------------------
            new_objective_value = self.f(new_x)
            if self.dimension == 1:
                self.regret += self.true_optimal - new_objective_value

            print ("iteration: {0}, new sample point: {1}, objective value: {2}, function value: {3}, average regret: {4}, beta: {5}".format(self.T, new_x, new_fun, new_objective_value, self.regret / (self.T+1), beta_t))
            # print ("iteration: {0}, optimal sample point: {1}, optimal objective value: {2}, function value: {3}".format(self.T, optimal_x, optimal_fun, optimal_objective_value))

            self.sample_points = np.concatenate((self.sample_points, np.reshape(new_x, (1, self.dimension))))
            self.sample_values = np.concatenate((self.sample_values, np.reshape(new_objective_value, (1,1))))
            self.T = self.T + 1

    def predict(self, x):
        gpr = GaussianProcessRegressor(kernel=self.kernel, optimizer=None, alpha=self.gp_alpha)
        gpr.fit(self.sample_points, self.sample_values)

        mean, std = gpr.predict(np.reshape(x, (len(x), self.dimension)), return_std=True)
        return mean, std

    def get_beta_t(self):
        if self.linear:
            beta_t = 2 * np.log(2 * self.T**2 * np.pi**2 / (3 * self.delta)) + 2 * self.dimension * np.log(self.T**2 * self.dimension * self.b * self.upper_bound * np.sqrt(np.log(4 * self.dimension * self.a / self.delta)))
        else:
            gamma_t = np.power(np.log(self.T), self.dimension + 1)
            beta_t = 2 * self.B + 300 * np.power(np.log(self.T / self.delta),3) * gamma_t
        return beta_t


    def GPUCB_objective_value(self, gpr, beta_t, x):
        mean, std = gpr.predict(np.reshape(x, (len(x), self.dimension)), return_std=True)
        mean = np.reshape(mean, (len(x)))
        return mean + np.sqrt(beta_t) * std

    def plot_current_prediction(self):
        gpr = self.gpr
        beta_t = self.get_beta_t()
        mean, std = gpr.predict(np.reshape(self.X_, (len(self.X_), self.dimension)), return_std=True)
        mean = np.reshape(mean, (len(self.X_)))
        plt.plot(self.X_, mean)
        plt.fill_between(np.reshape(self.X_, len(self.X_)), mean - np.sqrt(beta_t) * std, mean + np.sqrt(beta_t) * std, alpha=0.2, color='g')
        plt.show()



class DecomposedGPUCB: # TODO
    def __init__(self, decomposition, kernelList, dimension, upper_bound, constraints, delta=0.05, a=1, b=1, gp_alpha=0.01, initial_point=None, X_=None, discrete=False):
        self.f = decomposition.get_function_value
        self.decomposition = decomposition
        self.J = decomposition.J
        self.kernelList = kernelList
        self.T = 1
        self.dimension = dimension
        self.upper_bound = upper_bound
        self.constraints = constraints
        self.delta = delta
        self.a = a
        self.b = b
        self.gp_alpha = gp_alpha
        self.discrete = discrete
        if np.shape(X_):
            grid_size = 1000
            self.X_ = np.reshape(np.linspace(0, upper_bound, grid_size), (grid_size, 1))
        else:
            self.X_ = X_
        self.gpr_list = None
        self.regret = 0

        # ---------------------- initial sample -----------------------
        if not initial_point and dimension == 1:
            initial_point = np.random.rand() * upper_bound
        else:
            raise ("No specified initial point!")

        self.sample_points = np.reshape([initial_point], (1, dimension))
        self.sample_sub_values = np.reshape(self.decomposition.get_subfunction_values([initial_point]), (1, self.decomposition.J))
        self.sample_values = np.reshape([self.decomposition.get_function_value([initial_point])], (1,1))
        self.bds = [(0, upper_bound) for i in range(dimension)]
        if dimension == 1:
            self.true_optimal = np.max([self.decomposition.get_function_value([x]) for x in self.X_])
            print ("True optimal: {0}".format(self.true_optimal))
        else:
            self.true_optimal = None

    def run(self, iterations): # TODO
        for iteration in range(iterations):
            # beta_t = 2 * np.log(2 * self.T**2 * np.pi**2 / (3 * self.delta)) + 2 * self.dimension * np.log(self.T**2 * self.dimension * self.b * self.upper_bound * np.sqrt(np.log(4 * self.dimension * self.a / self.delta)))
            beta_t = self.get_beta_t()

            gpr_list = []
            for i in range(self.decomposition.J):
                gpr = GaussianProcessRegressor(kernel=self.kernelList[i], optimizer=None, alpha=self.gp_alpha)
                gpr.fit(self.sample_points, self.sample_sub_values[:,i])
                gpr_list.append(gpr)
            self.gpr_list = gpr_list

            fn = lambda x: -self.GPUCB_objective_value(gpr_list, beta_t, x)
            # --------------------- optimization choice -------------------
            if self.discrete:
                objective_values = [fn(x) for x in self.X_]
                optimal_index = np.argmin(objective_values)
                optimal_x = self.X_[optimal_index]
                optimal_fun = objective_values[optimal_index]
                optimal_objective_value = self.f(optimal_x)
                new_x = optimal_x
                new_fun = optimal_fun

            else:
                #initial_point = self.sample_points[np.argmax(self.sample_values)]
                initial_point = np.random.rand() * self.upper_bound

                res = scipy.optimize.minimize(fn, initial_point, bounds=self.bds, constraints=self.constraints)
                new_x = res.x
                new_fun = res.fun

            # ------------------- suboptimal choice -----------------------
            new_subfunction_values = self.decomposition.get_subfunction_values(new_x)
            new_objective_value = self.decomposition.g(new_subfunction_values)
            if self.dimension == 1:
                self.regret += self.true_optimal - new_objective_value

            print ("iteration: {0}, new sample point: {1}, objective value: {2}, function value: {3}, average regret: {4}, beta: {5}".format(self.T, new_x, new_fun, new_objective_value, self.regret / (self.T+1), beta_t))
            # print ("iteration: {0}, optimal sample point: {1}, optimal objective value: {2}, function value: {3}".format(self.T, optimal_x, optimal_fun, optimal_objective_value))

            self.sample_points = np.concatenate((self.sample_points, np.reshape(new_x, (1, self.dimension))))
            self.sample_values = np.concatenate((self.sample_values, np.reshape(new_objective_value, (1,1))))
            self.sample_sub_values = np.concatenate((self.sample_sub_values, np.reshape(new_subfunction_values, (1,self.J))))
            self.T = self.T + 1

    def predict(self, x):
        gpr = GaussianProcessRegressor(kernel=self.kernel, optimizer=None, alpha=self.gp_alpha)
        gpr.fit(self.sample_points, self.sample_values)

        mean, std = gpr.predict(np.reshape(x, (len(x), self.dimension)), return_std=True)
        return mean, std

    def get_beta_t(self):
        beta_t = 2 * np.log(2 * self.T**2 * np.pi**2 / (3 * self.delta)) + 2 * self.dimension * np.log(self.T**2 * self.dimension * self.b * self.upper_bound * np.sqrt(np.log(4 * self.dimension * self.a / self.delta)))
        return beta_t


    def GPUCB_objective_value(self, gpr_list, beta_t, x): # TODO currently only work for one x
        assert(len(x) == 1)
        overall_variance = np.zeros(len(x))
        individual_mean_list = np.zeros(self.decomposition.J)
        coefficient_list = self.decomposition.get_coefficients(x) # used for computing the variance bound
        for i in range(self.J):
            individual_mean, individual_std = gpr_list[i].predict(np.reshape(x, (len(x), self.dimension)), return_std=True)
            individual_mean = np.reshape(individual_mean, (len(x)))
            individual_mean_list[i] = individual_mean
            overall_variance += np.square(coefficient_list[i]) * np.square(individual_std)

        overall_mean = self.decomposition.g(individual_mean_list)
        overall_std = np.sqrt(overall_variance)

        return overall_mean + np.sqrt(beta_t) * overall_std

    def plot_current_prediction(self):
        gpr = self.gpr
        beta_t = self.get_beta_t()
        mean, std = gpr.predict(np.reshape(self.X_, (len(self.X_), self.dimension)), return_std=True)
        mean = np.reshape(mean, (len(x)))
        plt.plot(self.X_, mean)
        plt.fill_between(self.X_, mean - np.sqrt(beta_t) * std, mean + np.sqrt(beta_t) * std, alpha=0.2, color='g')
        plt.show()


def smoothify(target, X_):
    def f(x):
        index = np.argmax(x < X_)
        grid_size = len(X_)
        value = (index / float(grid_size) - x) * target[index-1] + (1 - index / float(grid_size) + x) * target[index]
        return value
    return f

def randomize(f_orig, gp_alpha):
    def f(x):
        return f_orig(x) + np.random.rand() * gp_alpha
    return f

def maxDerivative(target, grid_size, upper_bound=1):
    max_index = np.argmax(target)
    min_index = np.argmin(target)

    difference = target[max_index] - target[min_index]
    max_derivative = difference / (np.abs(max_index - min_index) * upper_bound / float(grid_size))
    return max_derivative

if __name__ == "__main__":
    upper_bound = 1
    grid_size = 1000
    kernel = Matern(length_scale=0.2)
    gp_alpha = 0.00001
    dimension = 1
    random_seed = np.random.randint(1000)
    #budget_constraint = 0.3
    #constraints = ({"type": "eq", "fun": lambda x: np.sum(x) - budget_constraint})
    constraints = None

    X_ = np.reshape(np.linspace(0, upper_bound, grid_size), (grid_size, 1))

    gp = GaussianProcessRegressor(kernel=kernel, optimizer=None, alpha=gp_alpha)
    target = gp.sample_y(X_, 1, random_seed)

    f = smoothify(target, X_)

    GPUCBsolver = GPUCB(f, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, X_=X_)
    #GPUCBsolver.run(1)




