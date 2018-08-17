import numpy as np
import pandas as pd

if __name__ == "__main__":
    directory = "result/"
    total_run = 200
    total_count = 5

    decomposed_min_list = []
    GPUCB_min_list = []
    decomposed_mean_list = []
    GPUCB_mean_list = []
    for i in range(30):
        try:
            decomposed_result = pd.read_csv(directory + 'decomposedGPUCB_result_0811_{0}.csv'.format(i), index_col=0).values / (total_run * total_count)
            GPUCB_result = pd.read_csv(directory + 'GPUCB_result_0811_{0}.csv'.format(i), index_col=0).values / (total_run * total_count)

            decomposed_result_min = np.min(decomposed_result)
            decomposed_result_mean = np.mean(decomposed_result)

            GPUCB_result_min = np.min(GPUCB_result)
            GPUCB_result_mean = np.mean(GPUCB_result)
            
            decomposed_min_list.append(decomposed_result_min)
            decomposed_mean_list.append(decomposed_result_mean)
            GPUCB_min_list.append(GPUCB_result_min)
            GPUCB_mean_list.append(GPUCB_result_mean)
        except:
            continue



