import numpy as np
import matplotlib.pyplot as plt

import argparse
import pickle
import pandas

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Decomposed GPUCB and GPUCB comparison')
    parser.add_argument('-n', '--name', required=True, help='Input the name to save')
    parser.add_argument('-s', '--scale_down', default=5, help='Input the scale down factor')
    parser.add_argument('-a', '--a', required=True, help='Input the a value')
    parser.add_argument('-b', '--b', required=True, help='Input the b value')
    parser.add_argument('-d', '--domain', required=True, help='Input the domain: flu, synthetic, or weather')


    args = parser.parse_args()
    filename = "{0}_scale{1}_a{2}_b{3}".format(args.name, args.scale_down, args.a, args.b)

    if args.domain == "synthetic":
        output_path = "synthetic/linear/"
    elif args.domain == "flu":
        # output_path = "flu/new_result/"
        output_path = "flu/result0825/"
    elif args.domain == "weather":
        output_path = "weather/new_result/"

    (GPUCB_regret_list, decomposedGPUCB_regret_list, EI_regret_list, decomposedEI_regret_list, POI_regret_list, decomposedPOI_regret_list) = pickle.load(open(output_path+"regret_list_{0}.p".format(filename), "rb"))
    
    total_run = GPUCB_regret_list.shape[1]

    average_GPUCB_regret_list = np.mean(GPUCB_regret_list, axis=0)
    average_decomposedGPUCB_regret_list = np.mean(decomposedGPUCB_regret_list, axis=0)
    average_EI_regret_list = np.mean(EI_regret_list, axis=0)
    average_decomposedEI_regret_list = np.mean(decomposedEI_regret_list, axis=0)
    average_POI_regret_list = np.mean(POI_regret_list, axis=0)
    average_decomposedPOI_regret_list = np.mean(decomposedPOI_regret_list, axis=0)

    average_GPUCB_regret = [np.mean(average_GPUCB_regret_list[:i+1]) for i in range(total_run)]
    average_decomposedGPUCB_regret = [np.mean(average_decomposedGPUCB_regret_list[:i+1]) for i in range(total_run)]
    average_EI_regret = [np.mean(average_EI_regret_list[:i+1]) for i in range(total_run)]
    average_decomposedEI_regret = [np.mean(average_decomposedEI_regret_list[:i+1]) for i in range(total_run)]
    average_POI_regret = [np.mean(average_POI_regret_list[:i+1]) for i in range(total_run)]
    average_decomposedPOI_regret = [np.mean(average_decomposedPOI_regret_list[:i+1]) for i in range(total_run)]

    #plt.yscale('log')
    plt.title(output_path)
    plt.plot(range(1, total_run+1), average_GPUCB_regret, 'b', label="GPUCB")
    plt.plot(range(1, total_run+1), average_decomposedGPUCB_regret, 'r', label="Decomposed")
    plt.plot(range(1, total_run+1), average_EI_regret, 'y', label="EI")
    plt.plot(range(1, total_run+1), average_POI_regret, 'g', label="POI")
    plt.legend()
    plt.xlabel('iterations')
    plt.ylabel('regret')
    #plt.show()
    plt.savefig(output_path + "visualize_{0}.png".format(filename))

    f = open(output_path+"summary_{0}.csv".format(filename), "w")
    f.write("GPUCB, " + ", ".join([str(x) for x in average_GPUCB_regret]) + "\n")
    f.write("decomposedGPUCB, " + ", ".join([str(x) for x in average_decomposedGPUCB_regret]) + "\n")
    f.write("EI, " + ", ".join([str(x) for x in average_EI_regret]) + "\n")
    f.write("decomposedEI, " + ", ".join([str(x) for x in average_decomposedEI_regret]) + "\n")
    f.write("POI, " + ", ".join([str(x) for x in average_POI_regret]) + "\n")
    f.write("decomposedPOI, " + ", ".join([str(x) for x in average_decomposedPOI_regret]) + "\n")
    f.close()

    print("improve ratio: {0}".format(average_GPUCB_regret[-1] / average_decomposedGPUCB_regret[-1]))
