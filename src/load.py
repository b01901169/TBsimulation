import numpy as np
import pickle

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

if __name__ == "__main__":
    file_index = "0124_null"
    dir_path = "/home/kai/Dropbox/USC/publication/TBsimulationCode/"

    fig = plt.figure()

    # =================== high fidelity model ==================
    print "============================ high fidelity ============================="

    # ------------------- ploting 3D surface ------------------
    high_infected, high_population = pickle.load(open(dir_path + "data/high_fidelity_{0}.data".format(file_index), "rb"))
    low_infected, low_population = pickle.load(open(dir_path + "data/low_fidelity_{0}.data".format(file_index), "rb"))

    high_infected = np.array(high_infected)
    high_population = np.array(high_population)
    high_percentage = high_infected / high_population

    x_length, y_length = high_infected.shape
    x_start_year = 30 # 30 years old started
    x_end_year = x_length - 50 # 60 years old
    y_year_length = y_length / 12
    start_year = 1996 - 130
    x1 = [i for i in range(x_start_year, x_end_year) for j in range(130, y_year_length)]
    y1 = [start_year + j for i in range(x_start_year, x_end_year) for j in range(130, y_year_length)]
    data_length = (x_end_year - x_start_year) * (y_year_length - 130)
    z1_infected = [high_infected[x1[i]][(y1[i] - start_year) * 12] for i in range(data_length)]
    z1_population = [high_population[x1[i]][(y1[i] - start_year) * 12] for i in range(data_length)]
    z1 = [high_percentage[x1[i]][(y1[i] - start_year) * 12] for i in range(data_length)]

    ax1 = fig.add_subplot(1, 3, 1, projection='3d')
    ax1.set_xlabel("age")
    ax1.set_ylabel("year")
    ax1.set_zlabel("# infected")

    ax1.plot_trisurf(x1, y1, z1, cmap=cm.coolwarm, linewidth=0.2, antialiased=True)
    #plt.show()


    #"""
    #"""
    # =================== low fidelity model ==================
    print "============================ low fidelity ============================="

    # ------------------- ploting 3D surface ------------------
    low_infected = np.array(low_infected)
    low_population = np.array(low_population)
    low_percentage = low_infected / low_population

    x_length, y_length = low_infected.shape
    start_age = 30
    start_year = 1995
    x2 = [i + start_age for i in range(x_length) for j in range(y_length)]
    y2 = [j + start_year for i in range(x_length) for j in range(y_length)]
    data_length = x_length * y_length
    low_infected = np.asarray(low_infected)
    z2_infected = [low_infected[x2[i] - start_age][y2[i] - start_year] for i in range(data_length)]
    z2_population = [low_population[x2[i] - start_age][y2[i] - start_year] for i in range(data_length)]
    z2 = [low_percentage[x2[i] - start_age][y2[i] - start_year] for i in range(data_length)]

    ax2 = fig.add_subplot(1, 3, 2, projection='3d')
    ax2.set_xlabel("age")
    ax2.set_ylabel("year")
    ax2.set_zlabel("# infected")

    ax2.plot_trisurf(x2, y2, z2, cmap=cm.coolwarm, linewidth=0.2, antialiased=True)

    # ========= difference between low and high fidelity ==========
    z3 = np.array(z2) - np.array(z1)

    ax3 = fig.add_subplot(1, 3, 3, projection='3d')
    ax3.set_xlabel("age")
    ax3.set_ylabel("year")
    ax3.set_zlabel("# infected")
    ax3.plot_trisurf(x1, y1, z3, cmap=cm.coolwarm, linewidth=0.2, antialiased=True)

    plt.show()

    total_variation = np.sum(np.abs(z3))

    print "file index: {0}".format(file_index)
    print "total variation: {0}".format(total_variation)
    print "simulation total infected {0}".format(np.sum(high_infected[30:60, range(130*12, 155*12, 12)]))
    print "approximated total infected {0}".format(np.sum(low_infected, 1))
    print "approximated total infected {0}".format(np.sum(low_infected))


    #"""

