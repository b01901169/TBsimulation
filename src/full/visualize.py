import numpy as np
import matplotlib.pyplot as plt

import argparse
import pickle
import pandas

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Decomposed GPUCB and GPUCB comparison')
    parser.add_argument('-n', '--name', help='Input the name to save')

    args = parser.parse_args()
    filename = args.name

    #output_path = "synthetic/linear/"
    output_path = "flu/result/"
    #output_path = "weather/result/"

    (GPUCB_regret_list, decomposed_regret_list, EI_regret_list, POI_regret_list) = pickle.load(open(output_path+"regret_list_{0}.p".format(filename), "rb"))
    
    GPUCB_scores = pandas.read_csv(output_path+"GPUCB_result_{0}.csv".format(filename), index_col=0)
    decomposedGPUCB_scores = pandas.read_csv(output_path+"decomposedGPUCB_result_{0}.csv".format(filename), index_col=0)
    EI_scores = pandas.read_csv(output_path+"EI_result_{0}.csv".format(filename), index_col=0)
    POI_scores = pandas.read_csv(output_path+"POI_result_{0}.csv".format(filename), index_col=0)

    GPUCB_best_index = np.unravel_index(GPUCB_scores.values.argmax(), GPUCB_scores.shape)
    decomposedGPUCB_best_index = np.unravel_index(decomposedGPUCB_scores.values.argmax(), decomposedGPUCB_scores.shape)
    EI_best_index = np.unravel_index(EI_scores.values.argmax(), EI_scores.shape)
    POI_best_index = np.unravel_index(POI_scores.values.argmax(), POI_scores.shape)

    total_run = GPUCB_regret_list.shape[3]

    average_GPUCB_regret_list = np.mean(GPUCB_regret_list, axis=2)
    average_decomposed_regret_list = np.mean(decomposed_regret_list, axis=2)
    average_EI_regret_list = np.mean(EI_regret_list, axis=2)
    average_POI_regret_list = np.mean(POI_regret_list, axis=2)

    average_GPUCB_regret = [np.mean(average_GPUCB_regret_list[GPUCB_best_index][:i+1]) for i in range(total_run)]
    average_decomposed_regret = [np.mean(average_decomposed_regret_list[decomposedGPUCB_best_index][:i+1]) for i in range(total_run)]
    average_EI_regret = [np.mean(average_EI_regret_list[EI_best_index][:i+1]) for i in range(total_run)]
    average_POI_regret = [np.mean(average_POI_regret_list[POI_best_index][:i+1]) for i in range(total_run)]

    #plt.yscale('log')
    plt.plot(range(1, total_run+1), average_GPUCB_regret, 'b')
    plt.plot(range(1, total_run+1), average_decomposed_regret, 'r')
    plt.plot(range(1, total_run+1), average_EI_regret, 'g')
    plt.plot(range(1, total_run+1), average_POI_regret, 'y')
    plt.show()
    #plt.savefig(output_path + "visualize_{0}.png".format(filename))

    f = open(output_path+"summary_{0}.csv".format(filename), "w")
    f.write("GPUCB, " + ", ".join([str(x) for x in GPUCB_regret_list[GPUCB_best_index]]) + "\n")
    f.write("decomposedGPUCB, " + ", ".join([str(x) for x in decomposed_regret_list[decomposedGPUCB_best_index]]) + "\n")
    f.close()

    print("improve ratio: {0}".format(average_GPUCB_regret[-1] / average_decomposed_regret[-1]))
