from CA import CA
from EM import EM
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time


def accuracy(actual_classification_vector, attempt_classification_vector):
    correct = 0
    for i in range(0, len(actual_classification_vector)):
        if actual_classification_vector[i] == attempt_classification_vector[i]:
            correct += 1
    return correct/len(actual_classification_vector)

def plot(InData, NumClust1, OutCenter1, NumClust2=None, OutCenter2=None, Dim=2, color=("red","black")):
    fig = plt.figure()
    ax = plt.axes()
    x = InData.iloc[:,0]
    y = InData.iloc[:,1]
    ax.scatter(x, y,s=1)
    
    centerCoord1 = np.zeros(Dim)
    for i in range(0, NumClust1):
        for j in range(0, Dim):
            centerCoord1[j] = OutCenter1[i][j]
        ax.scatter(centerCoord1[0], centerCoord1[1], color=color[0])
    
    if(NumClust2 and OutCenter2):
        centerCoord2 = np.zeros(Dim)
        for i in range(0, NumClust2):
            for j in range(0, Dim):
                centerCoord2[j] = OutCenter2[i][j]
            ax.scatter(centerCoord2[0], centerCoord2[1], color=color[1])
    plt.show()

if __name__ == '__main__':
    # Read data (InCenters.csv, InData.csv only things that are required)
    InCenters = pd.read_csv('data/training/random-centers-40.csv', header=None)
    InData = pd.read_csv('data/training/a1.csv', header=None)
    NumClust=25
    Iterations=50
    NumOfVectors=3000
    InCenters = InCenters[0:NumClust]
    actual_classification_file = pd.read_csv("data/training/a1-ga.csv", header=None)
    actual_classification_vector = actual_classification_file[0].tolist()

    # --------------------------CA-----------------------------------------
    print("Starting CA Algorithm:")
    CA_start_time = time.time_ns()
    CA_NumClust, CA_OutCenter, CA_Classifications = CA.CA(InData, InCenters, MaximumIt=Iterations, NumClust=NumClust, NumOfVectors=NumOfVectors)
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

    plot(InData, CA_NumClust, CA_OutCenter)