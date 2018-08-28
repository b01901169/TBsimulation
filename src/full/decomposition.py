import numpy as np
import scipy
from scipy.stats import norm

from matplotlib import pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import (RBF, Matern, RationalQuadratic,
                                              ExpSineSquared, DotProduct,
                                              ConstantKernel, WhiteKernel)

from sklearn.utils.validation import check_X_y, check_array

N_Y = True # Normalize Y
REPEAT_NUMBER = 3

class Decomposition:
    def __init__(self, J, fList, g, gList):
        # J: number of decomposed subfunctions
        # fList: list (size J) of subfunctions. The function is either a lambda function or a callable function with access to f(x)
        self.J = J
        self.fList = fList
        self.g = g
        self.gList = gList

        
    def get_function_value(self, x):
        # input: x (d dimension array)
        # output: function value f(x) = g(f1(x), f2(x), ..., fJ(x))
        # assert(type(x) == float or type(x) == int or len(x) == 1)

        subfunction_values = np.zeros(self.J)
        for i in range(self.J):
            subfunction_values[i] = self.fList[i](x)

        g_value = self.g(subfunction_values)
        return g_value


    def get_subfunction_values(self, x):
        # input: x (d dimension array)
        # output: J values representing the values of individual subfunctions
        # assert(type(x) == float or type(x) == int or len(x) == 1)

        subfunction_values = np.zeros(self.J)
        for i in range(self.J):
            subfunction_values[i] = self.fList[i](x)

        return subfunction_values

    def get_coefficients(self, x): # currently only allow for size 1
        # assert(type(x) == float or type(x) == int or len(x) == 1)
        coefficients = np.zeros((self.J, 1))
        for i in range(self.J):
            coefficients[i] = self.gList[i](x)

        return coefficients

class GPUCB:
    def __init__(self, f, kernel, dimension, upper_bound, constraints, delta=0.05, a=1, b=1, gp_alpha=0.01, initial_point=None, X_=None, discrete=False, linear=True, B=10, lower_bound=0, optimization_method=None, initial_point_generator=None, true_optimal=None, optimize_kernel=False, scale_down_factor=5, maxmin="max"):
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
        self.optimization_method = optimization_method
        self.initial_point_generator = initial_point_generator
        self.optimize_kernel = optimize_kernel
        self.scale_down_factor = scale_down_factor
        self.maxmin = maxmin
        if np.shape(X_):
            self.X_ = X_
            self.grid_size = len(X_)
        else:
            grid_size = 1000
            self.grid_size = grid_size
            self.X_ = np.reshape(np.linspace(0, upper_bound, grid_size), (grid_size, 1))
        self.gpr = None
        # ---------------------- regret records -----------------------
        self.regret = 0
        self.regret_list = []

        # ----------------------- true optimal ------------------------
        if true_optimal is not None:
            self.true_optimal = true_optimal
        elif discrete or dimension == 1:
            if self.maxmin == "max":
                self.true_optimal = np.max([f(x) for x in self.X_])
            elif self.maxmin == "min":
                self.true_optimal = np.min([f(x) for x in self.X_])
            print ("True optimal: {0}".format(self.true_optimal))
        else:
            self.true_optimal = None

        # ---------------------- initial sample -----------------------
        if initial_point is not None:
            initial_point = initial_point
        elif dimension == 1:
            initial_point = np.random.rand() * upper_bound
        elif discrete:
            initial_point = X_[np.random.randint(len(X_))]
        else:
            raise Exception("No specified initial point!")

        # if np.size(initial_point) == dimension:
        self.sample_points = np.reshape(initial_point, (1, dimension))
        self.sample_values = np.reshape([self.f(initial_point)], (1,1))
        self.bds = [(lower_bound, upper_bound) for i in range(dimension)]
        # else:
        #     self.sample_points = np.reshape(initial_point, (len(initial_point), dimension))
        #     self.sample_values = np.reshape([self.f(initial_point[i]) for i in range(len(initial_point))], (len(initial_point),1))
        #     self.bds = [(lower_bound, upper_bound) for i in range(dimension)]

    def run(self, iterations):
        for iteration in range(iterations):
            # beta_t = 2 * np.log(2 * self.T**2 * np.pi**2 / (3 * self.delta)) + 2 * self.dimension * np.log(self.T**2 * self.dimension * self.b * self.upper_bound * np.sqrt(np.log(4 * self.dimension * self.a / self.delta)))
            beta_t = self.get_beta_t()

            if self.optimize_kernel:
                gpr = GaussianProcessRegressor(kernel=self.kernel, alpha=self.gp_alpha, normalize_y=N_Y)
            else:
                gpr = GaussianProcessRegressor(kernel=self.kernel, optimizer=None, alpha=self.gp_alpha, normalize_y=N_Y)
            self.gpr = gpr
            gpr.fit(self.sample_points, self.sample_values)

            if self.maxmin == "max":
                fn = lambda x: -self.GPUCB_objective_value(gpr, beta_t, x, self.maxmin)
            elif self.maxmin == "min":
                fn = lambda x: self.GPUCB_objective_value(gpr, beta_t, x, self.maxmin)
            else:
                raise(Exception("Not implemented method!!"))
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
                new_x_list = []
                new_fun_list = []
                for j in range(REPEAT_NUMBER):
                    if self.initial_point_generator == "currentMax":
                        print("current max: {0}".format(np.max(self.sample_values)))
                        initial_point = self.sample_points[np.argmax(self.sample_values)]
                    elif self.initial_point_generator == "currentMin":
                        print("current min: {0}".format(np.min(self.sample_values)))
                        initial_point = self.sample_points[np.argmin(self.sample_values)]
                    elif self.initial_point_generator is not None:
                        initial_point = self.initial_point_generator()
                    else:
                        raise Exception("No initial point generator!!")

                    if self.optimization_method is None:
                        res = scipy.optimize.minimize(fn, initial_point, bounds=self.bds, constraints=self.constraints)
                    else:
                        res = scipy.optimize.minimize(fn, initial_point, bounds=self.bds, constraints=self.constraints, method=self.optimization_method)
                    # print (res)
                    new_x_list.append(res.x)
                    new_fun_list.append(res.fun)
                new_x = new_x_list[np.argmin(new_fun_list)]
                new_fun = new_fun_list[np.argmin(new_fun_list)]

            # ------------------- suboptimal choice -----------------------
            new_objective_value = self.f(new_x)
            if self.true_optimal != None:
                if self.maxmin == "max":
                    individual_regret = self.true_optimal - new_objective_value
                elif self.maxmin == "min":
                    individual_regret = - self.true_optimal + new_objective_value

                self.regret += individual_regret
                self.regret_list.append(individual_regret)

            format_new_x = ["{0:.2f}".format(x) for x in new_x]
            print ("iteration: {0}, new sample point: {1}, objective value: {2:.4f}, function value: {3:.4f}, average regret: {4:.4f}, beta: {5:.4f}".format(self.T, format_new_x, float(new_fun), new_objective_value, self.regret / (self.T), beta_t))
            # print ("iteration: {0}, optimal sample point: {1}, optimal objective value: {2}, function value: {3}".format(self.T, optimal_x, optimal_fun, optimal_objective_value))

            self.sample_points = np.concatenate((self.sample_points, np.reshape(new_x, (1, self.dimension))))
            self.sample_values = np.concatenate((self.sample_values, np.reshape(new_objective_value, (1,1))))
            self.T = self.T + 1

    def predict(self, x):
        gpr = GaussianProcessRegressor(kernel=self.kernel, optimizer=None, alpha=self.gp_alpha, normalize_y=N_Y)
        gpr.fit(self.sample_points, self.sample_values)

        mean, std = gpr.predict(np.reshape(x, (len(x), self.dimension)), return_std=True)
        return mean, std

    def get_beta_t(self): # TODO scale down by 5
        if self.discrete:
            beta_t = 2 * np.log(self.grid_size * (self.T ** 2) * (np.pi ** 2) / (6 * self.delta))
        elif self.linear:
            beta_t = 2 * np.log(2 * self.T**2 * np.pi**2 / (3 * self.delta)) + 2 * self.dimension * np.log(self.T**2 * self.dimension * self.b * self.upper_bound * np.sqrt(np.log(4 * self.dimension * self.a / self.delta)))
        else:
            gamma_t = np.power(np.log(self.T), self.dimension + 1)
            beta_t = 2 * self.B + 300 * np.power(np.log(self.T / self.delta),3) * gamma_t
        return beta_t / self.scale_down_factor


    def GPUCB_objective_value(self, gpr, beta_t, x, maxmin):
        x_len = int(x.size / self.dimension)
        if len(np.shape(x)) == 1:
            x = np.reshape(x, (x_len, self.dimension))
        try:
            check_array(x)
        except:
            return np.array([np.nan] * x_len)
            if maxmin == "max":
                return np.array([-np.inf]*x.shape[0])
            elif maxmin == "min":
                return np.array([np.inf]*x.shape[0])
        mean, std = gpr.predict(x, return_std=True)
        mean = np.reshape(mean, (x_len))
        if maxmin == "max":
            return mean + np.sqrt(beta_t) * std
        elif maxmin == "min":
            return mean - np.sqrt(beta_t) * std
        else:
            raise(Exception("Not implemented method!"))

    def plot_current_prediction(self):
        gpr = self.gpr
        beta_t = self.get_beta_t()
        mean, std = gpr.predict(np.reshape(self.X_, (len(self.X_), self.dimension)), return_std=True)
        mean = np.reshape(mean, (len(self.X_)))
        plt.plot(self.X_, mean)
        plt.fill_between(np.reshape(self.X_, len(self.X_)), mean - np.sqrt(beta_t) * std, mean + np.sqrt(beta_t) * std, alpha=0.2, color='g')
        plt.show()



class DecomposedGPUCB: # TODO
    def __init__(self, decomposition, kernelList, dimension, upper_bound, constraints, delta=0.05, a=1, b=1, gp_alpha=None, initial_point=None, X_=None, discrete=False, lower_bound=0, optimization_method=None, initial_point_generator=None, true_optimal=None, optimize_kernel=False, method="GPUCB", scale_down_factor=5, maxmin="max"):
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
        assert(len(gp_alpha) == self.J)
        self.gp_alpha = gp_alpha
        self.discrete = discrete
        self.optimization_method = optimization_method
        self.initial_point_generator = initial_point_generator
        self.optimize_kernel = optimize_kernel
        self.method = method
        self.scale_down_factor = scale_down_factor
        self.maxmin = maxmin
        if np.shape(X_):
            self.X_ = X_
            grid_size = X_.shape[0]
            self.grid_size = grid_size
        else:
            grid_size = 1000
            self.grid_size = grid_size
            self.X_ = np.reshape(np.linspace(0, upper_bound, grid_size), (grid_size, 1))
        self.gpr_list = None

        # ---------------------- regret records -----------------------
        self.regret = 0
        self.regret_list = []

        # ----------------------- true optimal ------------------------
        if true_optimal is not None:
            self.true_optimal = true_optimal
        elif discrete or (dimension == 1):
            if self.maxmin == "max":
                self.true_optimal = np.max([self.f(x) for x in self.X_])
            elif self.maxmin == "min":
                self.true_optimal = np.min([self.f(x) for x in self.X_])
            print ("True optimal: {0}".format(self.true_optimal))
        else:
            self.true_optimal = None
        # ---------------------- initial sample -----------------------
        if initial_point is not None:
            initial_point = initial_point
        elif dimension == 1:
            initial_point = np.random.rand() * upper_bound
        elif discrete:
            initial_point = X_[np.random.randint(len(X_))]
        else:
            raise Exception("No specified initial point!")

        if discrete:
            coefficient_list = np.zeros((self.J, grid_size))
            for i in range(grid_size):
                coefficient_list[:,i] = decomposition.get_coefficients(X_[i]).reshape(self.J)
            self.coefficient_list = coefficient_list
        else:
            coefficient_list = None
            self.coefficient_list = coefficient_list

        # if np.size(initial_point) == dimension:
        self.sample_points = np.reshape(initial_point, (1, dimension))
        self.sample_sub_values = np.reshape(self.decomposition.get_subfunction_values(initial_point), (1, self.decomposition.J))
        self.sample_values = np.reshape([self.f(initial_point)], (1,1))
        self.bds = [(lower_bound, upper_bound) for i in range(dimension)]
        # else:
        #     self.sample_points = np.reshape(initial_point, (len(initial_point), dimension))
        #     self.sample_sub_values = np.reshape([self.decomposition.get_subfunction_values(initial_point[i]) for i in range(len(initial_point))], (1, self.decomposition.J))
        #     self.sample_values = np.reshape([self.f(initial_point[i]) for i in range(len(initial_point))], (len(initial_point),1))
        #     self.bds = [(lower_bound, upper_bound) for i in range(dimension)]

    def run(self, iterations): # TODO
        for iteration in range(iterations):
            # beta_t = 2 * np.log(2 * self.T**2 * np.pi**2 / (3 * self.delta)) + 2 * self.dimension * np.log(self.T**2 * self.dimension * self.b * self.upper_bound * np.sqrt(np.log(4 * self.dimension * self.a / self.delta)))
            beta_t = self.get_beta_t()

            gpr_list = []
            for i in range(self.decomposition.J):
                if self.optimize_kernel:
                    gpr = GaussianProcessRegressor(kernel=self.kernelList[i], alpha=self.gp_alpha[i], normalize_y=N_Y)
                else:
                    gpr = GaussianProcessRegressor(kernel=self.kernelList[i], optimizer=None, alpha=self.gp_alpha[i], normalize_y=N_Y)
                gpr.fit(self.sample_points, self.sample_sub_values[:,i])
                gpr_list.append(gpr)
            self.gpr_list = gpr_list

            # --------------------- optimization choice -------------------
            if self.discrete:
                if self.maxmin == "max":
                    fn = lambda x: -self.GPUCB_objective_value(gpr_list, beta_t, x, maxmin=self.maxmin, whole_data=True)
                elif self.maxmin == "min":
                    fn = lambda x: self.GPUCB_objective_value(gpr_list, beta_t, x, maxmin=self.maxmin, whole_data=True)
                else:
                    raise(Exception("Not implemented method!!"))
                #fn = lambda x: -self.GPUCB_objective_value(gpr_list, beta_t, x)
                objective_values = fn(self.X_)
                optimal_index = np.argmin(objective_values)
                optimal_x = self.X_[optimal_index]
                optimal_fun = objective_values[optimal_index]
                optimal_objective_value = self.f(optimal_x)
                new_x = optimal_x
                new_fun = optimal_fun

            else:
                if self.maxmin == "max":
                    fn = lambda x: -self.GPUCB_objective_value(gpr_list, beta_t, x, self.maxmin)
                elif self.maxmin == "min":
                    fn = lambda x: self.GPUCB_objective_value(gpr_list, beta_t, x, self.maxmin)
                else:
                    raise(Exception("Not implemented method!!"))
                #initial_point = self.sample_points[np.argmax(self.sample_values)]
                new_x_list = []
                new_fun_list = []
                for j in range(REPEAT_NUMBER):
                    if self.initial_point_generator == "currentMax":
                        print("current max: {0}".format(np.max(self.sample_values)))
                        initial_point = self.sample_points[np.argmax(self.sample_values)]
                    elif self.initial_point_generator == "currentMin":
                        print("current min: {0}".format(np.min(self.sample_values)))
                        initial_point = self.sample_points[np.argmin(self.sample_values)]
                    elif self.initial_point_generator is not None:
                        initial_point = self.initial_point_generator()
                    else:
                        raise Exception("No initial point generator!!")

                    if self.optimization_method is None:
                        res = scipy.optimize.minimize(fn, initial_point, bounds=self.bds, constraints=self.constraints)
                    else:
                        res = scipy.optimize.minimize(fn, initial_point, bounds=self.bds, constraints=self.constraints, method=self.optimization_method)
                    new_x_list.append(res.x)
                    new_fun_list.append(res.fun)
                new_x = new_x_list[np.argmin(new_fun_list)]
                new_fun = new_fun_list[np.argmin(new_fun_list)]

            # ------------------- suboptimal choice -----------------------
            new_subfunction_values = self.decomposition.get_subfunction_values(new_x)
            new_objective_value = self.decomposition.g(new_subfunction_values)
            if self.true_optimal != None:
                if self.maxmin == "max":
                    individual_regret = self.true_optimal - new_objective_value
                elif self.maxmin == "min":
                    individual_regret = - self.true_optimal + new_objective_value

                self.regret += individual_regret
                self.regret_list.append(individual_regret)

            format_new_x = ["{0:.2f}".format(x) for x in new_x]
            print ("iteration: {0}, new sample point: {1}, objective value: {2:.4f}, function value: {3:.4f}, average regret: {4:.4f}, beta: {5:.4f}".format(self.T, format_new_x, float(new_fun), new_objective_value, self.regret / (self.T), beta_t))
            # print ("iteration: {0}, optimal sample point: {1}, optimal objective value: {2}, function value: {3}".format(self.T, optimal_x, optimal_fun, optimal_objective_value))

            self.sample_points = np.concatenate((self.sample_points, np.reshape(new_x, (1, self.dimension))))
            self.sample_values = np.concatenate((self.sample_values, np.reshape(new_objective_value, (1,1))))
            self.sample_sub_values = np.concatenate((self.sample_sub_values, np.reshape(new_subfunction_values, (1,self.J))))
            self.T = self.T + 1

    def predict(self, x): # TODO
        gpr = GaussianProcessRegressor(kernel=self.kernel, optimizer=None, alpha=self.gp_alpha, normalize_y=N_Y)
        gpr.fit(self.sample_points, self.sample_values)

        mean, std = gpr.predict(np.reshape(x, (len(x), self.dimension)), return_std=True)
        return mean, std

    def get_beta_t(self): # TODO scale down
        if self.discrete:
            beta_t = 2 * np.log(self.grid_size * (self.T ** 2) * (np.pi ** 2) / (6 * self.delta))
        else:
            beta_t = 2 * np.log(2 * self.T**2 * np.pi**2 / (3 * self.delta)) + 2 * self.dimension * np.log(self.T**2 * self.dimension * self.b * self.upper_bound * np.sqrt(np.log(4 * self.dimension * self.a / self.delta)))
        return beta_t / self.scale_down_factor


    def GPUCB_objective_value(self, gpr_list, beta_t, x, maxmin, whole_data=False, verbose=False): # TODO currently only work for one x
        if whole_data == False:
            assert(len(x) == self.dimension) # single data
            if len(np.shape(x)) == 1:
                x = np.reshape(x, (1, self.dimension))
            try:
                check_array(x)
            except:
                print('error:', x)
                if maxmin == "max":
                    return -np.inf
                elif maxmin == "min":
                    return np.inf

            overall_variance = 0
            individual_mean_list = np.zeros(self.J)
            coefficient_list = self.decomposition.get_coefficients(x) # used for computing the variance bound
            for i in range(self.J):
                individual_mean, individual_std = gpr_list[i].predict(x, return_std=True)
                individual_mean = np.reshape(individual_mean, (1))
                individual_mean_list[i] = individual_mean
                overall_variance += np.square(coefficient_list[i]) * np.square(individual_std)

            overall_mean = self.decomposition.g(individual_mean_list)
            overall_std = np.sqrt(overall_variance)
            if verbose:
                print("overall mean: {0}, overall std: {1}".format(overall_mean, overall_std))

            #return overall_mean + np.sqrt(beta_t) * overall_std

        else: # x.shape = (grid_size, dimension)
            assert(x.shape[1] == self.dimension)
            try:
                check_array(x)
            except:
                print('error:', x)
                if maxmin == "max":
                    return np.array([-np.inf]*x.shape[0])
                elif maxmin == "min":
                    return np.array([np.inf]*x.shape[0])

            grid_size = x.shape[0]
            overall_variance = np.zeros(grid_size)
            individual_mean_list = np.zeros((self.J, grid_size))
            coefficient_list = self.coefficient_list
            for i in range(self.J):
                individual_mean, individual_std = gpr_list[i].predict(x, return_std=True)
                # individual_mean, individual_std = gpr_list[i].predict(np.reshape(x, (1, self.dimension)), return_std=True)
                # individual_mean = np.reshape(individual_mean, (1))
                individual_mean_list[i] = individual_mean
                overall_variance += np.square(coefficient_list[i]) * np.square(individual_std)

            overall_mean = np.array([self.decomposition.g(individual_mean_list[:,i]) for i in range(grid_size)])
            overall_std = np.sqrt(overall_variance)

            #return overall_mean + np.sqrt(beta_t) * overall_std
        if maxmin == "max":
            best_f = np.max(self.sample_values)
            tmp_delta = overall_mean - best_f
        elif maxmin == "min":
            best_f = np.min(self.sample_values)
            tmp_delta = best_f - overall_mean

        if self.method == "GPUCB":
            if maxmin == "max":
                return overall_mean + np.sqrt(beta_t) * overall_std
            elif maxmin == "min":
                return overall_mean - np.sqrt(beta_t) * overall_std
            else:
                raise(Exception("Not implemented method!!"))

        elif self.method == "EI":
            x_len = int(x.size / self.dimension)
            tmp_delta_plus = np.max([tmp_delta, np.zeros(x_len)], axis=0)
            EI = tmp_delta_plus + overall_std * norm.pdf(tmp_delta/overall_std) - np.abs(tmp_delta) * norm.cdf(- np.abs(tmp_delta) / overall_std)
            assert(np.sum(EI < 0) == 0)
            return EI
        elif self.method == "POI":
            POI = norm.cdf(tmp_delta/overall_std)
            return POI
        else:
            raise(Exception("Not Implemented Method {0}!!".format(self.method)))

            # TODO

    def plot_current_prediction(self):
        gpr = self.gpr
        beta_t = self.get_beta_t()
        mean, std = gpr.predict(np.reshape(self.X_, (len(self.X_), self.dimension)), return_std=True)
        mean = np.reshape(mean, (len(x)))
        plt.plot(self.X_, mean)
        plt.fill_between(self.X_, mean - np.sqrt(beta_t) * std, mean + np.sqrt(beta_t) * std, alpha=0.2, color='g')
        plt.show()


class Improvement:
    def __init__(self, f, kernel, dimension, upper_bound, constraints, gp_alpha, method, X_=None, discrete=False, lower_bound=0, optimization_method=None, initial_point=None, initial_point_generator=None, true_optimal=None, optimize_kernel=False, maxmin="max"):
        self.f = f
        self.kernel = kernel
        self.T = 1
        self.dimension = dimension
        self.upper_bound = upper_bound
        self.constraints = constraints
        self.gp_alpha = gp_alpha
        self.initial_point = initial_point
        self.discrete = discrete
        self.method = method
        self.optimization_method = optimization_method
        self.initial_point_generator = initial_point_generator
        self.optimize_kernel = optimize_kernel
        self.maxmin = maxmin
        if np.shape(X_):
            self.X_ = X_
        else:
            grid_size = 1000
            self.X_ = np.reshape(np.linspace(0, upper_bound, grid_size), (grid_size, 1))

        # ---------------------- regret records -----------------------
        self.regret = 0
        self.regret_list = []

        # ----------------------- true optimal ------------------------
        if true_optimal is not None:
            self.true_optimal = true_optimal
        elif discrete or (dimension == 1):
            if self.maxmin == "max":
                self.true_optimal = np.max([f(x) for x in self.X_])
            elif self.maxmin == "min":
                self.true_optimal = np.min([f(x) for x in self.X_])
            print ("True optimal: {0}".format(self.true_optimal))
        else:
            self.true_optimal = None

        # ---------------------- initial sample -----------------------
        if initial_point is not None:
            initial_point = initial_point
        elif dimension == 1:
            initial_point = np.random.rand() * upper_bound
        elif discrete:
            initial_point = X_[np.random.randint(len(X_))]
        else:
            raise Exception("No specified initial point!")

        # if np.size(initial_point) == dimension:
        self.sample_points = np.reshape(initial_point, (1, dimension))
        self.sample_values = np.reshape([self.f(initial_point)], (1,1))
        self.bds = [(lower_bound, upper_bound) for i in range(dimension)]
        # else:
        #     self.sample_points = np.reshape(initial_point, (len(initial_point), dimension))
        #     self.sample_values = np.reshape([self.f(initial_point[i]) for i in range(len(initial_point))], (len(initial_point),1))
        #     self.bds = [(lower_bound, upper_bound) for i in range(dimension)]


    def run(self, iterations):
        for iteration in range(iterations):
            if self.optimize_kernel:
                gpr = GaussianProcessRegressor(kernel=self.kernel, alpha=self.gp_alpha, normalize_y=N_Y)
            else:
                gpr = GaussianProcessRegressor(kernel=self.kernel, optimizer=None, alpha=self.gp_alpha, normalize_y=N_Y)
            self.gpr = gpr
            gpr.fit(self.sample_points, self.sample_values)

            fn = lambda x: -self.expected_improvement_value(gpr, x)
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
                if self.initial_point_generator == "currentMax":
                    print("current max: {0}".format(np.max(self.sample_values)))
                    initial_point = self.sample_points[np.argmax(self.sample_values)]
                elif self.initial_point_generator == "currentMin":
                    print("current min: {0}".format(np.min(self.sample_values)))
                    initial_point = self.sample_points[np.argmin(self.sample_values)]
                elif self.initial_point_generator is not None:
                    initial_point = self.initial_point_generator()
                else:
                    raise Exception("No initial point generator!!")

                if self.optimization_method is None:
                    res = scipy.optimize.minimize(fn, initial_point, bounds=self.bds, constraints=self.constraints)
                else:
                    res = scipy.optimize.minimize(fn, initial_point, bounds=self.bds, constraints=self.constraints, method=self.optimization_method)
                # print (res)

                new_x = res.x
                new_fun = res.fun

            # ------------------- suboptimal choice -----------------------
            new_objective_value = self.f(new_x)
            if self.true_optimal != None:
                if self.maxmin == "max":
                    individual_regret = self.true_optimal - new_objective_value
                elif self.maxmin == "min":
                    individual_regret = - self.true_optimal + new_objective_value

                self.regret += individual_regret
                self.regret_list.append(individual_regret)

            format_new_x = ["{0:.2f}".format(x) for x in new_x]
            print ("iteration: {0}, new sample point: {1}, objective value: {2:.4f}, function value: {3:.4f}, average regret: {4:.4f}".format(self.T, format_new_x, float(new_fun), new_objective_value, self.regret / (self.T)))
            # print ("iteration: {0}, optimal sample point: {1}, optimal objective value: {2}, function value: {3}".format(self.T, optimal_x, optimal_fun, optimal_objective_value))

            self.sample_points = np.concatenate((self.sample_points, np.reshape(new_x, (1, self.dimension))))
            self.sample_values = np.concatenate((self.sample_values, np.reshape(new_objective_value, (1,1))))
            self.T = self.T + 1

    def expected_improvement_value(self, gpr, x):
        x_len = int(x.size / self.dimension)
        mean, std = gpr.predict(np.reshape(x, (x_len, self.dimension)), return_std=True)
        mean = np.reshape(mean, (x_len))
        if self.maxmin == "max":
            best_f = np.max(self.sample_values)
            delta = mean - best_f
        elif self.maxmin == "min":
            best_f = np.min(self.sample_values)
            delta = best_f - mean

        if self.method == "EI":
            delta_plus = np.max([delta, np.zeros(x_len)], axis=0)
            EI = delta_plus + std * norm.pdf(delta/std) - np.abs(delta) * norm.cdf(- np.abs(delta) / std)
            assert(np.sum(EI < 0) == 0)
            return EI
        elif self.method == "POI":
            POI = norm.cdf(delta/std)
            return POI
        else:
            raise(Exception("Not Implemented Method: {0}".format(self.method)))
        

def smoothify(target, X_):
    def f(x):
        index = np.argmax(x < X_)
        grid_size = len(X_)
        value = (index / float(grid_size) - x) * target[index-1] + (1 - index / float(grid_size) + x) * target[index]
        return value
    return f

def randomizify(f_orig, gp_alpha):
    def f(x):
        return f_orig(x) + (np.random.rand()-0.5) * gp_alpha
    return f

def maxDerivative(target, grid_size, upper_bound=1):
    max_index = np.argmax(target)
    min_index = np.argmin(target)

    difference = target[max_index] - target[min_index]
    max_derivative = difference / (np.linalg.norm(max_index - min_index) * upper_bound / float(grid_size))
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

    gp = GaussianProcessRegressor(kernel=kernel, optimizer=None, alpha=gp_alpha, normalize_y=N_Y)
    target = gp.sample_y(X_, 1, random_seed)

    f = smoothify(target, X_)

    GPUCBsolver = GPUCB(f, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, X_=X_)
    #GPUCBsolver.run(1)




