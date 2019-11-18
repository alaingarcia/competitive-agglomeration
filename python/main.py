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
"""
def plot(NumClust1, OutCenter1, NumClust2, OutCenter2, InData, Dim=2, color=("red","black")):
    fig = plt.figure()
    ax = plt.axes()
    x = InData.iloc[:,0]
    y = InData.iloc[:,1]
    ax.scatter(x, y)
    
    centerCoord1 = np.zeros(Dim)
    for i in range(0, NumClust1):
        for j in range(0, Dim):
            centerCoord1[j] = OutCenter1[i][j]
        ax.scatter(centerCoord1[0], centerCoord1[1], color=color[0])

    centerCoord2 = np.zeros(Dim)
    for i in range(0, NumClust2):
        for j in range(0, Dim):
            centerCoord2[j] = OutCenter2[i][j]
        ax.scatter(centerCoord2[0], centerCoord2[1], color=color[1])
    plt.show()
"""
if __name__ == '__main__':
    # Read data (InCenters.csv, InData.csv only things that are required)
    InCenters = pd.read_csv('data/InCenters.csv', header=None)
    InData = pd.read_csv('data/InData.csv', header=None)
    actual_classification_file = pd.read_csv("data/ActualClassification.csv", header=None)
    actual_classification_vector = actual_classification_file[0].tolist()

    # --------------------------CA-----------------------------------------
    CA_start_time = time.time_ns()
    CA_NumClust, CA_OutCenter, CA_Classifications = CA.CA(InData, InCenters, MaximumIt=50)
    CA_end_time = time.time_ns()

    # Print results
    print("CA final number of cluster: {}".format(CA_NumClust))
    print("CA cluster centers: {}".format(CA_OutCenter))
    print("CA classification vector: {}".format(CA_Classifications))

    # Convert ns to ms
    CA_time = (CA_end_time-CA_start_time)/1e6


    CA_accuracy = accuracy(actual_classification_vector,CA_Classifications)

    print("CA time: {} ms".format(CA_time))
    print("CA accuracy: {}%\n".format(CA_accuracy*100))

    # --------------------------EM-----------------------------------------

    EM_start_time = time.time_ns()
    EM_NumClust,EM_OutCenter, EM_Classifications  = EM.ExpectationMaximization(InData, InCenters, rep=50, pandas=True).info()
    EM_end_time = time.time_ns()

    print("EM final number of cluster: {}".format(EM_NumClust))
    print("EM cluster centers: {}".format(EM_OutCenter))
    print("EM classification vector: {}".format(EM_Classifications))

    EM_time = (EM_end_time-EM_start_time)/1e6

    EM_accuracy = accuracy(actual_classification_vector,EM_Classifications)

    print("EM time: {} ms".format(EM_time))
    print("EM accuracy: {}%\n\n".format(EM_accuracy*100))

    # PLOT
    #plot(CA_NumClust, CA_OutCenter, EM_NumClust, EM_OutCenter, InData)

