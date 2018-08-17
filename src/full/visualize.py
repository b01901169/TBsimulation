import numpy as np
import matplotlib.pyplot as plt

import argparse
import pickle

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Decomposed GPUCB and GPUCB comparison')
    parser.add_argument('-n', '--name', help='Input the name to save')

    args = parser.parse_args()
    filename = args.name

    input_path = "weather/result/"

    (GPUCB_regret_list, decomposed_regret_list) = pickle.load(open(input_path+"regret_list_{0}.p".format(filename), "r"))

