import numpy as np
import pandas as pd

if __name__ == "__main__":
    path = "./server0901/"
    f_output = open(path + "summary.csv", "w")
    for kernel_type in ["RBF", "RQ", "Matern"]:
        for J in [3,5,7,10]:
            for metric in ["rmse", "report"]:
                filename_0901 = path + "quantity_{0}_0901_new_{1}_J{2}.csv".format(metric, kernel_type, J)
                dataframe_0901 = pd.read_csv(filename_0901, header=None)
                # filename_0831 = path + "quantity_{0}_0831_{1}_J{2}.csv".format(metric, kernel_type, J)
                # dataframe_0831 = pd.read_csv(filename_0831, header=None)

                dataframe = np.concatenate([dataframe_0901.values[:,3:]], axis=0)[np.arange(0,150,3)+2]
                dataframe = np.array(dataframe, dtype=np.float32)
                dataframe = np.sqrt(dataframe) # originally mean square error, get the sqrt to get root-mean square error
                average_list = (np.mean(dataframe, axis=0) - 1) * 100
                median_list = (np.median(dataframe, axis=0) - 1) * 100

                f_output.write("{0}, {1}, {2},".format(metric, kernel_type, J) + ", ".join([str(x) for x in median_list]) + "\n")

    f_output.close()


