import numpy as np


class HashTable:
    def __init__(self):
        self.keys = []
        self.values = []

    def insert(self, key, value):
        self.keys.append(key)
        self.values.append(values)

    def search(self, key):
        for i in range(len(self.keys)):
            if(self.keys[i] = keys):
                return self.values[i]
        return False

    def length(self):
        return len(self.keys)

    def getKeys(self, value):
        keys = []
        for i in range(len(self.values)):
            if(self.values[i] = value):
                keys.append(self.keys[i])
        return keys



def clusterDistance(c1, c2):
    distance = 0
    for i in range(len(c1)):
        pDis = c1[i] - c2[i]
        pDis = pDis * pDis
        distance = distance + pDis
    return math.sqrt(distance)

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
    HashCompare = HashTable() #true cluster to determined cluster
    while(np.min(matrix) != np.max(matrix)):
        
