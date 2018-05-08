import numpy as np

from scipy.linalg import cholesky, cho_solve, solve_triangular

from sklearn.gaussian_process.kernels import RBF, ConstantKernel 
from sklearn.base import BaseEstimator, RegressorMixin, clone
from sklearn.utils import check_random_state
from sklearn.utils.validation import check_X_y, check_array
from sklearn.utils.deprecation import deprecated


class GaussianProcessRegressor:
    def __init__(self, kernel=None, alpha=1e-2,
                 normalize_y=False, copy_X_train=True, random_state=None):
        self.kernel = kernel
        self.alpha = alpha
        self.noise_level = alpha
        self.normalize_y = normalize_y
        self.copy_X_train = copy_X_train
        self.random_state = random_state 

    def fit_X(self, X):
        X = check_array(X)

        if self.kernel is None:  # Use an RBF kernel as default
            self.kernel_ = ConstantKernel(1.0, constant_value_bounds="fixed") \
                * RBF(1.0, length_scale_bounds="fixed")
        else:
            self.kernel_ = clone(self.kernel)


        self.X_train_ = np.copy(X) if self.copy_X_train else X

        # Precompute quantities required for predictions which are independent
        # of actual query points
        K = self.kernel_(self.X_train_)
        K[np.diag_indices_from(K)] += self.noise_level
        try:
            self.L_ = cholesky(K, lower=True)  # Line 2
        except np.linalg.LinAlgError as exc:
            exc.args = ("The kernel, %s, is not returning a "
                        "positive definite matrix. Try gradually "
                        "increasing the 'alpha' parameter of your "
                        "GaussianProcessRegressor estimator."
                        % self.kernel_,) + exc.args
            raise

        return self


    def fit_y(self, y):
        X, y = check_X_y(self.X_train_, y, multi_output=True, y_numeric=True)
        self.y_train_ = np.copy(y) if self.copy_X_train else y

        # Normalize target value
        if self.normalize_y:
            self._y_train_mean = np.mean(y, axis=0)
            # demean y
            y = y - self._y_train_mean
        else:
            self._y_train_mean = np.zeros(1)

        self.alpha_ = cho_solve((self.L_, True), self.y_train_)  # Line 3

        return self

    def predict(self, X, return_std=False, return_cov=False):
        
        if return_std and return_cov:
            raise RuntimeError(
                "Not returning standard deviation of predictions when "
                "returning full covariance.")

        X = check_array(X)

        K_trans = self.kernel_(X, self.X_train_)

        y_mean = K_trans.dot(self.alpha_)  # Line 4 (y_mean = f_star)
        y_mean = self._y_train_mean + y_mean  # undo normal.
        if return_cov:
            v = cho_solve((self.L_, True), K_trans.T)  # Line 5
            y_cov = self.kernel_(X) - K_trans.dot(v)  # Line 6
            return y_mean, y_cov
        elif return_std:
            # compute inverse K_inv of K based on its Cholesky
            # decomposition L and its inverse L_inv
            L_inv = solve_triangular(self.L_.T, np.eye(self.L_.shape[0]))
            K_inv = L_inv.dot(L_inv.T)
            # Compute variance of predictive distribution
            y_var = self.kernel_.diag(X)
            y_var -= np.einsum("ij,ij->i", np.dot(K_trans, K_inv), K_trans)

            # Check if any of the variances is negative because of
            # numerical issues. If yes: set the variance to 0.
            y_var_negative = y_var < 0
            if np.any(y_var_negative):
                warnings.warn("Predicted variances smaller than 0. "
                              "Setting those variances to 0.")
                y_var[y_var_negative] = 0.0
            return y_mean, np.sqrt(y_var)
        else:
            return y_mean

    def gradient(self, X, return_std=False): # only allow for one x
        X = check_array(X)

        K_trans = self.kernel_(X, self.X_train_)

        y_mean = K_trans.dot(self.alpha_)  # Line 4 (y_mean = f_star)
        y_mean = self._y_train_mean + y_mean

        alpha_X_train = self.alpha_.reshape(-1,1) * self.X_train_
        y_gradient = y_mean * X - K_trans.dot(alpha_X_train)

        if return_std: # TODO
            return y_gradient
        else:
            return y_gradient
        

    """
    def tensor_gradient(self, X, return_std=False):
        X = check_array(X)

        K_trans = self.kernel_(X, self.X_train_)

        y_mean = K_trans.dot(self.alpha_)  # Line 4 (y_mean = f_star)
        y_mean = self._y_train_mean + y_mean

        tensor_K_trans = 

        y_gradient = y_mean * X - K_trans.dot(alpha_X_train)

        if return_std: # TODO
            return y_gradient
        else:
            return y_gradient
    #"""
