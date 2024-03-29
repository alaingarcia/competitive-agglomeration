from EM import EM
import pandas as pd
import numpy as np
import time
from accuracy import matrixAccuracy

def EM_RUN(NumClust=25):
    #-------------------------small db-------------------------------------
    # InCenters = pd.read_csv('data/A1-Dataset/random-centers-40.csv', header=None)
    # InCenters = InCenters.to_numpy().tolist()
    #
    # InData_origin = pd.read_csv('data/A1-Dataset/a1.csv', header=None)
    # InData = InData_origin.to_numpy().tolist()
    #
    # Iterations=50
    # NumOfVectors=3000
    # InCenters = InCenters[0:NumClust]
    #
    # actual_classification_file = pd.read_csv("data/A1-Dataset/a1-ga.csv", header=None)
    # actual_classification_vector = actual_classification_file[0].tolist()
    #
    # actual_clusters = pd.read_csv("data/A1-Dataset/a1-ga-cb.csv", header=None)
    # actual_clusters = actual_clusters.to_numpy().tolist()
    #
    # trainIndexes = pd.read_csv("data/A1-Dataset/training_index.csv", header=None)
    # trainIndexes = trainIndexes[0].tolist()
    #--------------------------Large DB-------------------------------------

    InData = pd.read_csv('data/mnist_csv.csv', header=None)

    actual_classV = InData.iloc[:,-1]
    actual_classification_vector = actual_classV.tolist()

    tFile = open("./data/trainingMNIST.csv")
    trainIndexes = [int(i) for i in tFile.readlines()]
    tFile.close()

    InData = InData.iloc[:,:-1]
    InData = InData.to_numpy().tolist()
    InCenters = pd.read_csv('data/mnistCenterPoints.csv', header=None)
    InCenters = InCenters.to_numpy().tolist()

    actual_clusters = pd.read_csv("./data/averageTrue.csv")
    actual_clusters = actual_clusters.iloc[:,:-1]
    actual_clusters = actual_clusters.to_numpy().tolist()


    #--------------------------EM-----------------------------------------

    EM_train_start_time = time.time_ns()
    EM_train_NumClust,EM_train_OutCenter, EM_train_Classifications  = EM.ExpectationMaximization(InData, InCenters, rep=50, testSet=trainIndexes, classification=actual_classification_vector).info()
    EM_train_end_time = time.time_ns()
    EM_train_time = (EM_train_end_time-EM_train_start_time)/1e6
    EM_train_accuracy = matrixAccuracy(actual_clusters, actual_classification_vector, EM_train_OutCenter, EM_train_Classifications)

    EM_start_time = time.time_ns()
    EM_NumClust,EM_OutCenter, EM_Classifications  = EM.ExpectationMaximization(InData, InCenters, rep=50).info()
    EM_end_time = time.time_ns()
    EM_time = (EM_end_time-EM_start_time)/1e6
    EM_accuracy = matrixAccuracy(actual_clusters, actual_classification_vector, EM_OutCenter, EM_Classifications)

    out = open('EM_train_OutCenter.csv', 'w')
    for row in EM_train_OutCenter:
        for column in row:
            out.write('%d,' % column)
        out.write('\n')
    out.close()

    out = open('EM_OutCenter.csv', 'w')
    for row in EM_OutCenter:
        for column in row:
            out.write('%d,' % column)
        out.write('\n')
    out.close()

    out = open('EM_train_Classifications.csv', 'w')
    for row in EM_train_Classifications:
        out.write('%d,' % row)
        out.write('\n')
    out.close()

    out = open('EM_Classifications.csv', 'w')
    for row in EM_Classifications:
        out.write('%d,' % row)
        out.write('\n')
    out.close()

    #
    # # print("EM_trained cluster centers: {}".format(EM_train_OutCenter))
    # # print("EM_trained classification vector: {}".format(EM_train_Classifications))
    print("EM_trained final number of cluster: {}".format(EM_train_NumClust))
    print("EM_trained time: {} ms".format(EM_train_time))
    print("EM_trained accuracy: {}%\n\n".format(EM_train_accuracy*100))

    # print("EM cluster centers: {}".format(EM_OutCenter))
    # print("EM classification vector: {}".format(EM_Classifications))
    print("EM final number of cluster: {}".format(EM_NumClust))
    print("EM time: {} ms".format(EM_time))
    print("EM accuracy: {}%\n\n".format(EM_accuracy*100))


    #---------------determine matrix for assigning clusters--------------------
    return (EM_NumClust, EM_OutCenter)

if __name__ == '__main__':
    EM_RUN()
