import numpy as np
import math


class HashTable:
    def __init__(self):
        self.keys = []
        self.values = []

    def insert(self, key, value):
        self.keys.append(key)
        self.values.append(value)

    def search(self, key):
        for i in range(len(self.keys)):
            if(self.keys[i] == key):
                return self.values[i]
        return None

    def has(self, key):
        for i in range(len(self.keys)):
            if(self.keys[i] == key):
                return True
        return False

    def length(self):
        return len(self.keys)

    def getKeys(self, value):
        keys = []
        for i in range(len(self.values)):
            if(self.values[i] == value):
                keys.append(self.keys[i])
        if(keys == []):
            return None
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
    # 3.Match true cluster to determined cluster
    # 4.replace indexes from Determined Cluster Values to True Cluster values () make a mock classification list
    # 5.Determine Accuracy

    #1 - For each Determined cluster, find the index of the points that are attached to it
    indexCollector = []
    for i in range(len(attempt_clusters)): # i is the cluster
        for j in range(len(attempt_classV)): # j is the index of the points
            if(i == attempt_classV[j]): # if the classification of this point is the current cluster
                indexCollector.append(j)

    #2 - Creates and fills matrix
    matrix = np.zeros((len(attempt_clusters), len(actual_clusters)))
    for i in range(len(attempt_clusters)):
        for j in range(len(actual_clusters)):
            matrix[i][j] = clusterDistance(attempt_clusters[i],actual_clusters[j])

    #3 - Use matrix to check which clusters belongs to which other
    HashCompare = HashTable() # determined cluster to true cluster
    while(np.min(matrix) != np.max(matrix)):
        min_i = -1
        min_j = -1
        value = None
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                if(value == None or matrix[i][j] < value):
                    value = matrix[i][j]
                    min_i = i
                    min_j = j

        HashCompare.insert(min_i,min_j)
        newRow = [np.max(matrix)] * len(matrix[min_i])
        newCol = [np.max(matrix)] * len(matrix[:, min_j])
        matrix[min_i] = newRow
        matrix[:, min_j] = newCol

    #HashCompare would now have DETERMINED -> TRUE

    #4 - Create a new translated classification list to compare from filtered to actual_classV
    filtered_classification = [-1] * len(actual_classV) #make a list of [-1]
    for i in range(len(attempt_classV)):
        if(HashCompare.has(attempt_classV[i])):
            filtered_classification[i] = HashCompare.search(attempt_classV[i])

    #5 - Get accuracy for translated classifications
    correct = 0
    for i in range(len(actual_classV)):
        if actual_classV[i] == filtered_classification[i]:
            correct += 1
    return correct/len(actual_classV)
