from CA import CA
from EM import EM
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
import time


def accuracy(actual_classification_vector, attempt_classification_vector):
    correct = 0
    for i in range(0, len(actual_classification_vector)):
        if actual_classification_vector[i] == attempt_classification_vector[i]:
            correct += 1
    return correct/len(actual_classification_vector)


if __name__ == '__main__':
    # Read data (InCenters.csv, InData.csv only things that are required)
    #Centers = pd.read_csv('data/', header=None)
    InData = pd.read_csv('data/mnist_csv.csv', header=None)
    classV = InData.iloc[:,-1]
    InData = InData.iloc[:,:-1]

    # actual_classification_file = pd.read_csv("data/ActualClassification.csv", header=None)
    # actual_classification_vector = actual_classification_file[0].tolist()

    # --------------------------CA-----------------------------------------

    # EM_start_time = time.time_ns()
    # EM_NumClust,EM_OutCenter, EM_Classifications  = EM.ExpectationMaximization(InData, InCenters, rep=50, pandas=True).info()
    # EM_end_time = time.time_ns()
    #
    # print("EM final number of cluster: {}".format(EM_NumClust))
    # print("EM cluster centers: {}".format(EM_OutCenter))
    # print("EM classification vector: {}".format(EM_Classifications))
    #
    # EM_time = (EM_end_time-EM_start_time)/1e6
    #
    # EM_accuracy = accuracy(actual_classification_vector,EM_Classifications)
    #
    # print("EM time: {} ms".format(EM_time))
    # print("EM accuracy: {}%\n\n".format(EM_accuracy*100))
