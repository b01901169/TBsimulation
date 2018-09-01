import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
from sklearn.gaussian_process.kernels import (RBF, Matern, RationalQuadratic,
                                              ExpSineSquared, DotProduct,
                                              ConstantKernel, WhiteKernel,
                                              RationalQuadratic)
import argparse
import pickle

from decomposition import *


def g(function_values): # T_kelvin: temperature (kelvin), R: relative humidity, v_meter_sec: wind speed meter/second
    assert(len(function_values) == 3)
    T_kelvin = function_values[0]
    R = function_values[1]
    v_meter_sec = function_values[2]
    T = (T_kelvin - 273.15) * (1.8) + 32

    if T >= 70: # heat index case
        c1 = -42.38
        c2 = 2.049
        c3 = 10.14
        c4 = -0.2248
        c5 = -0.006838
        c6 = -0.05482
        c7 = 0.001228
        c8 = 0.0008528
        c9 = -0.00000199
        
        heat_index = c1 + c2 * T + c3 * R + c4 * T * R + c5 * T * T + c6 * R * R + c7 * T * T * R + c8 * T * R * R + c9 * T * T * R * R
        return heat_index

    elif T <= 50:
        v = v_meter_sec * 2.23684
        c1 = 35.74
        c2 = 0.6215
        c3 = -35.75
        c4 = 0.4275
        wind_chill = c1 + c2 * T + c3 * np.power(v, 0.16) + c4 * T * np.power(v, 0.16)

        return wind_chill

    else:
        return T

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Decomposed GPUCB and GPUCB comparison')
    parser.add_argument('-n', '--name', required=True, help='Input the name to save')
    parser.add_argument('-s', '--scale_down', default=1, help='Input the scale down factor')
    parser.add_argument('-iteration', '--iteration', default=300, help='Input the total iterations')
    parser.add_argument('-count', '--count', default=10, help='Input the total count')
    parser.add_argument('-a', '--a', required=True, help='Input the a value')
    parser.add_argument('-b', '--b', required=True, help='Input the b value')

    args = parser.parse_args()
    scale_down_factor = float(args.scale_down)
    filename = "{0}_scale{1}_a{2}_b{3}".format(args.name, args.scale_down, args.a, args.b)
    total_run = int(args.iteration)
    total_count = int(args.count)

    data_path = "weather/weather_14.json"
    output_path = "weather/new_result/"

    city_list = []
    timeframe_list = []
    temperature_list = []
    humidity_list = []
    wind_list = []
    perceived_list = []
    x_list = []
    y_list = []
    x_shift = 130
    y_shift = -20
    X_original = []
    for line in open(data_path, 'r'):
        weather_tmp = json.loads(line)
        if weather_tmp['city']['country'] == "US" and weather_tmp['city']['coord']['lon'] >= -130:
            # print(weather_tmp['city']['name'])
            city_list.append(weather_tmp['city']['name'])
            timeframe_list.append(weather_tmp['time'])

            T_kelvin = weather_tmp['main']['temp']
            humidity = weather_tmp['main']['humidity']
            wind_speed = weather_tmp['wind']['speed']

            x_coord = weather_tmp['city']['coord']['lon'] + x_shift
            y_coord = weather_tmp['city']['coord']['lat'] + y_shift

            temperature_list.append(T_kelvin)
            humidity_list.append(humidity)
            wind_list.append(wind_speed)
            perceived_list.append(g([T_kelvin, humidity, wind_speed]))
            x_list.append(x_coord)
            y_list.append(y_coord)

            X_original.append((x_coord, y_coord))

    X_original = np.array(X_original)
    x_list = np.array(x_list)
    y_list = np.array(y_list)

    humidity_list = np.array(humidity_list)
    temperature_list = np.array(temperature_list)
    wind_list = np.array(wind_list)

    # ============================= visualization ===========================
    # plt.figure(1)
    # plt.subplot(311)
    # sc1 = plt.scatter(x_list, y_list, c=humidity_list, cmap="Greens", marker='o')
    # plt.colorbar(sc1)
    # plt.subplot(312)
    # sc2 = plt.scatter(x_list, y_list, c=temperature_list, cmap="Reds", marker='o')
    # plt.colorbar(sc2)
    # plt.subplot(313)
    # sc3 = plt.scatter(x_list, y_list, c=perceived_list, cmap="hot", marker='o')
    # plt.colorbar(sc3)
    # 
    # plt.show()

    # ======================== individual function ==========================
    x2index = {}
    for i in range(len(X_original)):
        x = tuple(X_original[i])
        x2index[x] = i

    def temperatureFunction(x):
        i = x2index[tuple(x)]
        return temperature_list[i]

    def humidityFunction(x):
        i = x2index[tuple(x)]
        return humidity_list[i]

    def windFunction(x):
        i = x2index[tuple(x)]
        return wind_list[i]

    # ========================= experimental setting ========================
    J = 3
    dimension = 2
    upper_bound = 65
    constraints = None
    gp_alpha = 0.05
    gp_alpha_list = [0.04, 0.03, 0.005] # TODO gp alpha list
    delta = 0.05
    linear = True
    discrete = True
    optimize_kernel = False
    grid_size_original = len(temperature_list)

    # ============================= decomposition ===========================
    fList = [randomizify(temperatureFunction, gp_alpha_list[0]), randomizify(humidityFunction, gp_alpha_list[1]), randomizify(windFunction, gp_alpha_list[2])]

    function_bounds = np.ones(J)
    # for i in range(J):
    #     function_bounds[i] = np.mean(np.prod([np.abs(targetList[j]) for j in np.delete(np.arange(J), i)], axis=0))
    #gList = [lambda x: 0.5, lambda x: 0.15, lambda x: 0.35]
    gList = [lambda x: 1, lambda x: 0.3 if x[0] > 294.261 else 0, lambda x: 0 if x[0] > 283.15 else 0.3]

    decomposition = Decomposition(J, fList, g, gList)
    max_derivative_list = [maxDerivative(temperature_list, grid_size_original), maxDerivative(humidity_list, grid_size_original), maxDerivative(wind_list, grid_size_original)]

    # =================================== kernels determination =========================

    # kernelList = [13.8**2 * RBF(length_scale=13.5, length_scale_bounds=(2e-2, 2e2))     + 2.32**2 * Matern(length_scale=0.663, length_scale_bounds=(2e-2, 2e2)),
    #               17.0**2 * RBF(length_scale=4.35, length_scale_bounds=(2e-2, 2e2))     + 16.8**2 * Matern(length_scale=0.371, length_scale_bounds=(2e-2, 2e2)),
    #               1.23**2 * RBF(length_scale=8.00, length_scale_bounds=(2e-2, 2e2))     + 1.13**2 * Matern(length_scale=0.287, length_scale_bounds=(2e-2, 2e2))]
    # kernel =      54.6**2 * RBF(length_scale=32.5, length_scale_bounds=(2e-2, 2e2))     + 7.92**2 * Matern(length_scale=1.780, length_scale_bounds=(2e-2, 2e2))

    # kernelList = [11.6**2 * RBF(length_scale=18.10, length_scale_bounds=(2e-2, 2e2)),
    #               19.8**2 * RBF(length_scale=6.64, length_scale_bounds=(2e-2, 2e2)),
    #               2.44**2 * RBF(length_scale=6.01, length_scale_bounds=(2e-2, 2e2))]
    # kernel =      50.0**2 * RBF(length_scale=20.3, length_scale_bounds=(2e-2, 2e2))

    # kernelList = np.array(kernelList) #* 0.33
    # kernel = kernel #* 0.33

    # def initial_point_generator():
    #     initial_point = X_[np.random.randint(len(X_))]
    #     return initial_point

    #"""
    # ========================== experimental design ======================================
    a = float(args.a) 
    b = float(args.b) # * np.mean(max_derivative_list)
    print("a: {0}, b: {1}".format(a,b))
    maxmin = "max"

    # ================================ experimental records ===============================
    GPUCB_scores = np.zeros(total_count)
    decomposedGPUCB_scores = np.zeros(total_count)
    EI_scores = np.zeros(total_count)
    POI_scores = np.zeros(total_count)
    decomposedEI_scores = np.zeros(total_count)
    decomposedPOI_scores = np.zeros(total_count)

    GPUCB_regret_list = np.zeros((total_count, total_run))
    decomposedGPUCB_regret_list = np.zeros((total_count, total_run))
    EI_regret_list = np.zeros((total_count, total_run))
    POI_regret_list = np.zeros((total_count, total_run))
    decomposedEI_regret_list = np.zeros((total_count, total_run))
    decomposedPOI_regret_list = np.zeros((total_count, total_run))

    f_output = open(output_path + "scores_report_{0}.csv".format(filename), 'w')
    f_output.write("round, GPUCB, decomposed GPUCB, EI, decomposed EI, POI, decomposed POI")
    for count in range(total_count):
        # GPUCB_scores = np.zeros((a_count, b_count))
        # decomposedGPUCB_scores = np.zeros((a_count, b_count))
        # =================================== kernels fitting ===============================
        kernelList = []

        kernel_sample_size = 1000
        # kernel_sample_size = 1500
        tmp_x_list = np.zeros((kernel_sample_size, dimension))
        subfunction_values = np.zeros((kernel_sample_size, J))
        function_values = np.zeros((kernel_sample_size))
        #x_choices = np.random.binomial(1, size=grid_size, p=float(kernel_sample_size)/grid_size)
        x_choices = np.random.choice(grid_size_original, replace=False, size=kernel_sample_size)
        x_remaining = list(set(range(grid_size_original)) - set(x_choices))
        xList = X_original[x_choices]
        X_ = X_original[x_remaining]
        grid_size = len(X_)

        kernel_sample_size = sum(x_choices)
        print("grid size:", X_.shape)
        for i in range(len(xList)):
            x = xList[i]
            tmp_x_list[i] = x
            subfunction_values[i] = decomposition.get_subfunction_values(x)
            function_values[i] = decomposition.get_function_value(x)
        # ------------------------- whole function -----------------------------
        #gpr = GaussianProcessRegressor(kernel=1.0*Matern(length_scale=1, length_scale_bounds=(1, 2e2)), normalize_y=True)
        #gpr = GaussianProcessRegressor(kernel=1.0*RBF(length_scale=10, length_scale_bounds=(5, 20)) + 1.0*Matern(length_scale=1, length_scale_bounds=(2e-2, 5))
        #        + 1.0*RationalQuadratic(alpha=0.1, length_scale=1, length_scale_bounds=(2e-2, 1)), normalize_y=True)
        gpr = GaussianProcessRegressor(kernel=1.0*RBF(length_scale=10, length_scale_bounds=(10,10)) + 
                                              1.0*Matern(length_scale=1, length_scale_bounds=(2e-2, 5)) +
                                              1.0*RationalQuadratic(alpha=0.1, length_scale=1, length_scale_bounds=(2e-2, 1)), normalize_y=True)
        gpr.fit(tmp_x_list, function_values)
        kernel = gpr.kernel_

        # --------------------------- sub function -----------------------------
        for i in range(J):
            #gpr = GaussianProcessRegressor(kernel=1.0*Matern(length_scale=1, length_scale_bounds=(1, 2e2)), normalize_y=True)
            if i == 0:
                gpr = GaussianProcessRegressor(kernel=1.0*RBF(length_scale=10, length_scale_bounds=(10,10)) + 
                                                      1.0*Matern(length_scale=1, length_scale_bounds=(2e-2, 5), nu=1.5) +
                                                      1.0*RationalQuadratic(alpha=0.1, length_scale=1, length_scale_bounds=(2e-2, 1)), normalize_y=True)
            else:
                gpr = GaussianProcessRegressor(kernel=1.0*RBF(length_scale=10, length_scale_bounds=(10,10)) +
                                                      1.0*Matern(length_scale=1, length_scale_bounds=(2e-2, 5), nu=1.5) +
                                                      1.0*RationalQuadratic(alpha=0.1, length_scale=1, length_scale_bounds=(2e-2, 1)), normalize_y=True)
            gpr.fit(tmp_x_list, subfunction_values[:,i])
            kernelList.append(gpr.kernel_)
        #"""

        print("whole kernel: {0}".format(kernel))
        print("kernel list: {0}".format(kernelList))
        # =====================================================================================

        initial_point = X_[np.random.randint(grid_size)]
        f_output.write("\n{0}, ".format(count))

        print ("\nGPUCB count:{0}...".format(count))
        GPUCBsolver = GPUCB(decomposition.get_function_value, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, a=a, b=b, X_=X_, initial_point=initial_point, delta=delta, discrete=discrete, linear=linear, optimize_kernel=optimize_kernel, scale_down_factor=scale_down_factor, maxmin=maxmin) # linear arg only changes the beta_t used in exploration
        GPUCBsolver.run(total_run)
        GPUCB_scores[count] = GPUCBsolver.regret
        GPUCB_regret_list[count] = np.array(GPUCBsolver.regret_list)
        f_output.write("{0}, ".format(GPUCBsolver.regret))
        print("kernel: {0}".format(GPUCBsolver.gpr.kernel_))

        print ("\ndecomposed GPUCB count:{0}...".format(count))
        decomposedGPUCBsolver = DecomposedGPUCB(decomposition, kernelList, dimension, upper_bound, constraints, gp_alpha=gp_alpha_list, a=a, b=b, X_=X_, initial_point=initial_point, delta=delta, discrete=discrete, optimize_kernel=optimize_kernel, scale_down_factor=scale_down_factor, maxmin=maxmin)
        decomposedGPUCBsolver.run(total_run)
        decomposedGPUCB_scores[count] = decomposedGPUCBsolver.regret
        decomposedGPUCB_regret_list[count] = np.array(decomposedGPUCBsolver.regret_list)
        f_output.write("{0}, ".format(decomposedGPUCBsolver.regret))
        print("kernel: {0}".format([decomposedGPUCBsolver.gpr_list[i].kernel_ for i in range(J)]))

        print ("\nExpected Improvement count:{0}...".format(count))
        EIsolver = Improvement(decomposition.get_function_value, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, method="EI", X_=X_, initial_point=initial_point, discrete=discrete, optimize_kernel=optimize_kernel, maxmin=maxmin)
        EIsolver.run(total_run)
        EI_scores[count] = EIsolver.regret
        EI_regret_list[count] = np.array(EIsolver.regret_list)
        f_output.write("{0}, ".format(EIsolver.regret))

        # print ("\ndecomposed EI count:{0}...".format(count))
        # decomposedEIsolver = DecomposedGPUCB(decomposition, kernelList, dimension, upper_bound, constraints, gp_alpha=gp_alpha_list, a=a, b=b, X_=X_, initial_point=initial_point, delta=delta, discrete=discrete, optimize_kernel=optimize_kernel, method="EI", maxmin=maxmin)
        # decomposedEIsolver.run(total_run)
        # decomposedEI_scores[count] = decomposedEIsolver.regret
        # decomposedEI_regret_list[count] = np.array(decomposedEIsolver.regret_list)
        # f_output.write("{0}, ".format(decomposedEIsolver.regret))

        print ("\nProbability of Improvement count:{0}...".format(count))
        POIsolver = Improvement(decomposition.get_function_value, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, method="POI", X_=X_, initial_point=initial_point, discrete=discrete, optimize_kernel=optimize_kernel, maxmin=maxmin)
        POIsolver.run(total_run)
        POI_scores[count] = POIsolver.regret
        POI_regret_list[count] = np.array(POIsolver.regret_list)
        f_output.write("{0}, ".format(POIsolver.regret))

        # print ("\ndecomposed POI count:{0}...".format(count))
        # decomposedPOIsolver = DecomposedGPUCB(decomposition, kernelList, dimension, upper_bound, constraints, gp_alpha=gp_alpha_list, a=a, b=b, X_=X_, initial_point=initial_point, delta=delta, discrete=discrete, optimize_kernel=optimize_kernel, method="POI")
        # decomposedPOIsolver.run(total_run)
        # decomposedPOI_scores[count] = decomposedPOIsolver.regret
        # decomposedPOI_regret_list[count] = np.array(decomposedPOIsolver.regret_list)
        # f_output.write("{0}, ".format(decomposedPOIsolver.regret))

    f_output.close()
    pickle.dump((GPUCB_regret_list, decomposedGPUCB_regret_list, EI_regret_list, decomposedEI_regret_list, POI_regret_list, decomposedPOI_regret_list), open(output_path+"regret_list_{0}.p".format(filename), 'wb'))
