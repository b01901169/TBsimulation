import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
from sklearn.gaussian_process.kernels import (RBF, Matern, RationalQuadratic,
                                              ExpSineSquared, DotProduct,
                                              ConstantKernel, WhiteKernel)
import argparse
import pandas
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
    parser.add_argument('-n', '--name', help='Input the name to save')

    args = parser.parse_args()
    filename = args.name

    data_path = "weather/weather_14.json"
    output_path = "weather/result/"

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
    grid_size = len(temperature_list)

    # ============================= decomposition ===========================
    fList = [randomizify(temperatureFunction, gp_alpha_list[0]), randomizify(humidityFunction, gp_alpha_list[1]), randomizify(windFunction, gp_alpha_list[2])]
    kernelList = [RBF(length_scale=5), 0.5*RBF(length_scale=5), 0.1*RBF(length_scale=5)]
    kernel = 0.5*RBF(length_scale=5)

    function_bounds = np.ones(J)
    # for i in range(J):
    #     function_bounds[i] = np.mean(np.prod([np.abs(targetList[j]) for j in np.delete(np.arange(J), i)], axis=0))
    #gList = [lambda x: 0.5, lambda x: 0.15, lambda x: 0.35]
    gList = [lambda x: 1, lambda x: 0.3, lambda x: 0.7]

    decomposition = Decomposition(J, fList, g, gList, kernelList)
    max_derivative_list = [maxDerivative(temperature_list, grid_size), maxDerivative(humidity_list, grid_size), maxDerivative(wind_list, grid_size)]

    # ========================== experimental design ========================
    total_count = 10
    total_run = 300
    a_count = 5
    #a_list = np.array([0.05]) * np.mean(max_derivative_list)
    a_list = np.array([1e-5, 2e-5, 5e-5, 0.0001, 0.0002]) * np.mean(max_derivative_list)
    b_count = 5
    #b_list = np.array(np.arange(0.05, 0.06, 0.01))
    b_list = np.array(np.arange(0.01, 0.06, 0.01))

    GPUCB_scores = np.zeros((a_count, b_count))
    decomposedGPUCB_scores = np.zeros((a_count, b_count))
    GPUCB_regret_list = np.zeros((a_count, b_count, total_run))
    decomposed_regret_list = np.zeros((a_count, b_count, total_run))

    for count in range(total_count):
        # GPUCB_scores = np.zeros((a_count, b_count))
        # decomposedGPUCB_scores = np.zeros((a_count, b_count))

        for a_index in range(a_count):
            a = a_list[a_index]
            for b_index in range(b_count):
                b = b_list[b_index]

                initial_point = X_[np.random.randint(grid_size)]

                print ("\nGPUCB...")
                GPUCBsolver = GPUCB(decomposition.get_function_value, kernel, dimension, upper_bound, constraints, gp_alpha=gp_alpha, a=a, b=b, X_=X_, initial_point=initial_point, delta=delta, discrete=discrete, linear=linear) # linear arg only changes the beta_t used in exploration
                GPUCBsolver.run(total_run)
                GPUCB_scores[a_index, b_index] += GPUCBsolver.regret
                GPUCB_regret_list[a_index, b_index] += np.array(GPUCBsolver.regret_list)

                print ("\ndecomposed GPUCB")
                decomposedGPUCBsolver = DecomposedGPUCB(decomposition, kernelList, dimension, upper_bound, constraints, gp_alpha=gp_alpha_list, a=a, b=b, X_=X_, initial_point=initial_point, delta=delta, discrete=discrete)
                decomposedGPUCBsolver.run(total_run)
                decomposedGPUCB_scores[a_index, b_index] += decomposedGPUCBsolver.regret
                decomposed_regret_list[a_index, b_index] += np.array(decomposedGPUCBsolver.regret_list)

    GPUCB_df = pandas.DataFrame(data=GPUCB_scores, columns=b_list, index=a_list)
    decomposedGPUCB_df = pandas.DataFrame(data=decomposedGPUCB_scores, columns=b_list, index=a_list)

    GPUCB_df.to_csv(path_or_buf=output_path+'GPUCB_result_{0}.csv'.format(filename))
    decomposedGPUCB_df.to_csv(path_or_buf=output_path+'decomposedGPUCB_result_{0}.csv'.format(filename))

    decomposed_regret_list /= total_count
    GPUCB_regret_list /= total_count
    pickle.dump((GPUCB_regret_list, decomposed_regret_list), open(output_path+"regret_list_{0}.p".format(filename), 'wb'))
