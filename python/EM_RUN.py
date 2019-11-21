from CA import CA
from EM import EM
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
import time

def clusterDistance(c1, c2):
    distance = 0
    for i in range(len(c1)):
        pDis = c1[i] - c2[i]
        pDis = pDis * pDis
        distance = distance + pDis
    return math.sqrt(distance)


def accuracy(actual_classification_vector, attempt_classification_vector):
    correct = 0
    for i in range(0, len(actual_classification_vector)):
        if actual_classification_vector[i] == attempt_classification_vector[i]:
            correct += 1
    return correct/len(actual_classification_vector)

def matrixAccuracy(actual_clusters, actual_classV, attempt_clusters, attempt_classV):
    #Checklist:
    # 1.Get indexes for each Determined Cluster in personal code
    # 2.Create matrix to match a True Cluster to a Determined Cluster
    # 3.replace indexes from Determined Cluster Values to True Cluster values
    # 4.Determine Accuracy

    #1
    #For each Determined cluster, find the index of the points that are attached to it
    indexCollector = []
    for i in range(EM_NumClust): # i is the cluster
        for j in range(len(EM_Classification)): # j is the index of the points
            if(i == EM_Classification[j]): # if the classification of this point is the current cluster
                indexCollector.append(j)

    #2- Created and fills matrix
    matrix = np.zeros(len(actual_clusters), len(attempt_clusters))
    for i in range(len(actual_clusters)):
        for j in range(len(attempt_clusters)):
            matrix[i][j] = clusterDistance(actual_clusters[i],attempt_clusters[j])

    #3 - Use matrix to check which clusters belongs to which other
    matrixMax = np.max(matrix)





if __name__ == '__main__':
    # Read data (InCenters.csv, InData.csv only things that are required)
    InData = pd.read_csv('data/mnist_csv.csv', header=None)

    classV = InData.iloc[:,-1]
    actual_classification_vector = classV.tolist()

    InData = InData.iloc[:,:-1]

    InCenters = pd.read_csv('data/mnistCenterPoints.csv', header=None)

    tFile = open("./data/trainingMNIST.txt")
    trainIndexes = [int(i) for i in tFile.readlines()]
    tFile.close()

    #--------------------------EM-----------------------------------------

    # EM_train_start_time = time.time_ns()
    # EM_train_NumClust,EM_train_OutCenter, EM_train_Classifications  = EM.ExpectationMaximization(InData, InCenters, rep=50, testSet=trainIndexes, classification=actual_classification_vector, pandas=True).info()
    # EM_train_end_time = time.time_ns()
    # EM_train_time = (EM_train_end_time-EM_train_start_time)/1e6
    # EM_train_accuracy = accuracy(actual_classification_vector,EM_train_Classifications)

    # EM_start_time = time.time_ns()
    # EM_NumClust,EM_OutCenter, EM_Classifications  = EM.ExpectationMaximization(InData, InCenters, rep=50, pandas=True).info()
    # EM_end_time = time.time_ns()
    # EM_time = (EM_end_time-EM_start_time)/1e6
    # EM_accuracy = accuracy(actual_classification_vector,EM_Classifications)
    #
    # print("EM_trained final number of cluster: {}".format(EM_train_NumClust))
    # # print("EM_trained cluster centers: {}".format(EM_train_OutCenter))
    # # print("EM_trained classification vector: {}".format(EM_train_Classifications))
    #
    # print("EM_trained time: {} ms".format(EM_train_time))
    # print("EM_trained accuracy: {}%\n\n".format(EM_train_accuracy*100))
    #
    # print("EM final number of cluster: {}".format(EM_NumClust))
    # # print("EM cluster centers: {}".format(EM_OutCenter))
    # # print("EM classification vector: {}".format(EM_Classifications))
    #
    # print("EM time: {} ms".format(EM_time))
    # print("EM accuracy: {}%\n\n".format(EM_accuracy*100))


    #---------------determine matrix for assigning clusters--------------------
