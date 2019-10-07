import CA
from sys import argv
import pandas as pd
import numpy as np
from sklearn import datasets
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

MaximumIt = 10
EITA_0 = 5.000000
TAU = 10.000000
EXPINC = 20.000000
EXPDEC = 80.000000
MIN_SIMILARITY = 0.1
EPS = 1.0e-6
MAX_CENT_DIFF = 0.001

if __name__ == '__main__':
    '''NumClust = argv[1]
    NumOfVectors = argv[2]
    Dim = argv[3]
    Eita_0 = argv[4]
    MaximumIt = argv[5]
    TAU = argv[6]
    EXPINC = argv[7]
    EXPDEC = argv[8]
    m_p = argv[9] # power for distance
    ca_path = argv[10]'''

    NumClust = 4
    NumOfVectors = 50
    Dim = 3
    m_p = 3
 
    print("Start Reading")
    input = datasets.load_iris()

    # Initialize
    FeatVect = []
    for i in range(0, NumOfVectors):
        FeatVect.insert(i,CA.Feature_Info(np.zeros(NumClust), 
                        np.zeros(NumClust), np.zeros(Dim), None))
    # Populate
    for i in range(0, Dim):
        for j in range(0, NumOfVectors):
            FeatVect[j].dimen[i] = input.data[j,i]

    # Centers start at 0, 0
    # centers = np.zeros((NumClust, Dim))
    centers = np.array([[0.,0.,0.],[1.,1.,1.],[2.,2.,2.],[3.,3.,3.]])

    # MinPts?
    MinPts = np.zeros(MaximumIt)
    MinPts2 = np.zeros(MaximumIt)

    NumClust = CA.CA(FeatVect, NumOfVectors, Dim, NumClust, centers,
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
    ax = plt.axes(projection='3d')
    x = input.data[:,0]
    y = input.data[:,1]
    z = input.data[:,2]
    ax.scatter3D(x, y, z, c=input.target, cmap='Set1')
    
    centerCoord = np.zeros(Dim)
    for i in range(0, NumClust):
        for j in range(0, Dim):
            centerCoord[j] = centers[i][j]
        ax.scatter3D(centerCoord[0], centerCoord[1], centerCoord[2], color='black')
    plt.show()