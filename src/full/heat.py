import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
from sklearn.gaussian_process.kernels import (RBF, Matern, RationalQuadratic,
                                              ExpSineSquared, DotProduct,
                                              ConstantKernel, WhiteKernel)
import argparse
import pickle

from decomposition import *


def g(function_values): # T_kelvin: temperature (kelvin), R: relative humidity, v_meter_sec: wind speed meter/second
    assert(len(function_values) == 3)
    T_kelvin = function_values[0]
    R = function_values[1]
    v_meter_sec = function_values[2]
    T = (T_kelvin - 273.15) * (1.8) + 32

    if T >= 80: # heat index case
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
    elif T <= 50: # wind chill
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
    parser.add_argument('-s', '--scale_down', required=True, help='Input the scale down factor')
    parser.add_argument('-iteration', '--iteration', default=300, help='Input the total iterations')
    parser.add_argument('-count', '--count', default=10, help='Input the total count')

    args = parser.parse_args()
    scale_down_factor = float(args.scale_down)
    filename = "{0}_{1}".format(args.name, args.scale_down)
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
    X_ = []
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

            X_.append((x_coord, y_coord))

    X_ = np.array(X_)
    x_list = np.array(x_list)
    y_list = np.array(y_list)

    humidity_list = np.array(humidity_list)
    temperature_list = np.array(temperature_list)
    wind_list = np.array(wind_list)

    # ============================= visualization ===========================
    # plt.figure(1)
    # plt.subplot(311)
    # plt.scatter(x_list, y_list, c=humidity_list, cmap="Greens", marker='o')
    # plt.subplot(312)
    # plt.scatter(x_list, y_list, c=temperature_list, cmap="Reds", marker='o')
    # plt.subplot(313)
    # plt.scatter(x_list, y_list, c=perceived_list, cmap="hot", marker='o')
    # 
    # plt.show()

    # ======================== individual function ==========================
    x2index = {}
    for i in range(len(X_)):
        x = tuple(X_[i])
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
    gp_alpha = 1
    gp_alpha_list = [1, 10, 0.5] # TODO gp alpha list
    delta = 0.05
    linear = True
    discrete = True
    optimize_kernel = True
    grid_size = len(temperature_list)

    # ============================= decomposition ===========================
    fList = [randomizify(temperatureFunction, gp_alpha_list[0]), randomizify(humidityFunction, gp_alpha_list[1]), randomizify(windFunction, gp_alpha_list[2])]
    kernelList = [RBF(length_scale=5, length_scale_bounds=(1, 20)),
              0.5*RBF(length_scale=5, length_scale_bounds=(1, 20)),
              0.1*RBF(length_scale=5, length_scale_bounds=(1, 20))]
    kernel =  RBF(length_scale=5, length_scale_bounds=(1, 20))

    function_bounds = np.ones(J)
    # for i in range(J):
    #     function_bounds[i] = np.mean(np.prod([np.abs(targetList[j]) for j in np.delete(np.arange(J), i)], axis=0))
    #gList = [lambda x: 0.5, lambda x: 0.15, lambda x: 0.35]
    gList = [lambda x: 1, lambda x: 0.3, lambda x: 0.7]

    decomposition = Decomposition(J, fList, g, gList, kernelList)
    max_derivative_list = [maxDerivative(temperature_list, grid_size), maxDerivative(humidity_list, grid_size), maxDerivative(wind_list, grid_size)]

    # ========================== experimental design ========================
    a = np.mean(max_derivative_list)
    b = 0.5

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

        initial_point = X_[np.random.randint(grid_size)]
        f_output.write("\n{0}, ".format(count))

        print ("\nGPUCB count:{0}...".format(count))
        GPUCBsolver = GPUCB(decomposition.get_function_value, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, a=a, b=b, X_=X_, initial_point=initial_point, delta=delta, discrete=discrete, linear=linear, optimize_kernel=optimize_kernel, scale_down_factor=scale_down_factor) # linear arg only changes the beta_t used in exploration
        GPUCBsolver.run(total_run)
        GPUCB_scores[count] = GPUCBsolver.regret
        GPUCB_regret_list[count] = np.array(GPUCBsolver.regret_list)
        f_output.write("{0}, ".format(GPUCBsolver.regret))

        print ("\ndecomposed GPUCB count:{0}...".format(count))
        decomposedGPUCBsolver = DecomposedGPUCB(decomposition, kernelList, dimension, upper_bound, constraints, gp_alpha=gp_alpha_list, a=a, b=b, X_=X_, initial_point=initial_point, delta=delta, discrete=discrete, optimize_kernel=optimize_kernel, scale_down_factor=scale_down_factor)
        decomposedGPUCBsolver.run(total_run)
        decomposedGPUCB_scores[count] = decomposedGPUCBsolver.regret
        decomposedGPUCB_regret_list[count] = np.array(decomposedGPUCBsolver.regret_list)
        f_output.write("{0}, ".format(decomposedGPUCBsolver.regret))

        print ("\nExpected Improvement count:{0}...".format(count))
        EIsolver = Improvement(decomposition.get_function_value, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, method="EI", X_=X_, initial_point=initial_point, discrete=discrete, optimize_kernel=optimize_kernel)
        EIsolver.run(total_run)
        EI_scores[count] = EIsolver.regret
        EI_regret_list[count] = np.array(EIsolver.regret_list)
        f_output.write("{0}, ".format(EIsolver.regret))

        print ("\ndecomposed EI count:{0}...".format(count))
        decomposedEIsolver = DecomposedGPUCB(decomposition, kernelList, dimension, upper_bound, constraints, gp_alpha=gp_alpha_list, a=a, b=b, X_=X_, initial_point=initial_point, delta=delta, discrete=discrete, optimize_kernel=optimize_kernel, method="EI")
        decomposedEIsolver.run(total_run)
        decomposedEI_scores[count] = decomposedEIsolver.regret
        decomposedEI_regret_list[count] = np.array(decomposedEIsolver.regret_list)
        f_output.write("{0}, ".format(decomposedEIsolver.regret))

        print ("\nProbability of Improvement count:{0}...".format(count))
        POIsolver = Improvement(decomposition.get_function_value, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, method="POI", X_=X_, initial_point=initial_point, discrete=discrete, optimize_kernel=optimize_kernel)
        POIsolver.run(total_run)
        POI_scores[count] = POIsolver.regret
        POI_regret_list[count] = np.array(POIsolver.regret_list)
        f_output.write("{0}, ".format(POIsolver.regret))

        print ("\ndecomposed POI count:{0}...".format(count))
        decomposedPOIsolver = DecomposedGPUCB(decomposition, kernelList, dimension, upper_bound, constraints, gp_alpha=gp_alpha_list, a=a, b=b, X_=X_, initial_point=initial_point, delta=delta, discrete=discrete, optimize_kernel=optimize_kernel, method="POI")
        decomposedPOIsolver.run(total_run)
        decomposedPOI_scores[count] = decomposedPOIsolver.regret
        decomposedPOI_regret_list[count] = np.array(decomposedPOIsolver.regret_list)
        f_output.write("{0}, ".format(decomposedPOIsolver.regret))

    f_output.close()
    pickle.dump((GPUCB_regret_list, decomposedGPUCB_regret_list, EI_regret_list, decomposedEI_regret_list, POI_regret_list, decomposedPOI_regret_list), open(output_path+"regret_list_{0}.p".format(filename), 'wb'))
