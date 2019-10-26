import CA
from sys import argv
import pandas as pd
import numpy as np
from sklearn import datasets
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt



if __name__ == '__main__':
    MaximumIt = 50
    EITA_0 = 2.5
    TAU = 10
    EXPINC = 25
    EXPDEC = 35
    MIN_SIMILARITY = 0.1
    EPS = 1.0e-6
    MAX_CENT_DIFF = 0.001
    NumClust = 10
    NumOfVectors = 60
    Dim = 2
    m_p = 2
 
    print("Start Reading")
    InCenters = pd.read_csv('InCenters.csv', header=None)
    InData = pd.read_csv('InData.csv', header=None)
    MinClSize = np.zeros((50,2))

    # Initialize Data (from InData)
    FeatVect = []
    for i in range(0, NumOfVectors):
        FeatVect.insert(i,CA.Feature_Info(np.zeros(NumClust), 
                        np.zeros(NumClust), np.zeros(Dim), None))
    # Populate
    for i in range(0, Dim):
        for j in range(0, NumOfVectors):
            FeatVect[j].dimen[i] = InData.iloc[j,i]

    # Initialize Centers (random data from InCenters)
    centers = InCenters.values

    # MinPts are zero vectors ?
    MinPts = np.zeros(MaximumIt)
    MinPts2 = np.ones(MaximumIt)
    MinPts2[int(MaximumIt/2):MaximumIt] *= 5

    NumClust, Center = CA.CA(FeatVect, NumOfVectors, Dim, NumClust, centers,
                    MaximumIt, EITA_0, MinPts, MinPts2, TAU,
                    EXPINC, EXPDEC, m_p)

    print(NumClust, file=open('NumClust.txt', 'w'))

    clusterOut = open('OutCluster.txt', 'w')
    for i in range(0, NumOfVectors):
        clusterOut.write(str(FeatVect[i].cluster) + '\n')
    clusterOut.close()
    
    centerOut = open('OutCenters.txt', 'w')
    for i in range(0, NumClust):
        for j in range(0, Dim):
            centerOut.write(str(centers[i][j]) + ' ')
        centerOut.write('\n')
    centerOut.close()

    fig = plt.figure()
    ax = plt.axes()
    x = InData.iloc[:,0]
    y = InData.iloc[:,1]
    ax.scatter(x, y)
    
    centerCoord = np.zeros(Dim)
    for i in range(0, NumClust):
        for j in range(0, Dim):
            centerCoord[j] = centers[i][j]
        ax.scatter(centerCoord[0], centerCoord[1], color='black')
    plt.show()