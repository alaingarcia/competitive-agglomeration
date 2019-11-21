import math
import pandas as pd

class Point:
    def __init__(self, dimArray):
        self.features = dimArray
        self.cluster = -1

    def assignCluster(self, clusterIndex):
        self.cluster = clusterIndex

class Cluster:
    def __init__(self, mu):
        self.mu = mu

    def pointDistance(self, pointV): # Use euclidean distance to determine which cluster each point belongs to
        distance = 0
        for i in range(len(pointV.features)):
            pDis = pointV.features[i] - self.mu[i]
            pDis = pDis * pDis
            distance = distance + pDis
        return math.sqrt(distance)

    def __str__(self): #string version of a cluster
        return "Cluster contains mu values of: {0}".format(self.mu)

# holds points and means for each training set
class TestSet:
    def __init__(self, group, points, mu):
        self.group = group
        self.mu = mu
        self.points = points
        self.cluster = -1

    def clusterDistance(self, cluster):
        distance = 0
        for i in range(len(cluster.mu)):
            pDis = cluster.mu[i] - self.mu[i]
            pDis = pDis * pDis
            distance = distance + pDis
        return math.sqrt(distance)

class ExpectationMaximization:

    #Complete algorithm in init
    def __init__(self, data, clusterData, rep=100, testSet=None, classification=None, pandas=False):

        if(pandas):
            data, clusterData = self.pandasConvert(data, clusterData)

        if(len(data[0]) != len(clusterData[0])):
            print("ERROR: points and clusters do not have the same length")
            exit(0)

        if(testSet and classification): #testSet and classification has value
            self.trainingSets = self.trainingSetsCreation(data, testSet, classification)
        self.dataPoints = self.pointCreation(data)
        self.clusters = self.clusterCreation(clusterData)

        if(testSet and classification and len(self.clusters) < len(self.trainingSets)):
            print("ERROR: There are more groups than there are clusters")
            exit(0)

        for i in range(rep):
            self.E_step(self.clusters, self.dataPoints)
            self.M_step(self.clusters, self.dataPoints)
            if(testSet and classification):
                self.Train(self.clusters, self.dataPoints, self.trainingSets)
                self.M_step(self.clusters, self.dataPoints)
            print(i)

    #Expectation Step
    def E_step(self, clusters, data):
        # assign cluster to points
        for point in data:
            minValue = math.inf
            bestClus = -1
            for i in range(len(clusters)):
                test = clusters[i].pointDistance(point)
                if test < minValue:
                    minValue = test
                    bestClus = i
            point.assignCluster(bestClus)

        # count the number of points to each cluster
        pCount = [0] * len(clusters)
        for point in data:
            if(point.cluster != -1):
                pCount[point.cluster] = pCount[point.cluster] + 1
            else:
                    print("ERROR: Point missing cluster!!")
                    exit()

        # check if any clusters have less than 5% of available points
        clusterIndex = []
        redo = False
        for i in range(len(pCount)):
            if(pCount[i] < len(data) * (1/len(self.clusters)) * (4/5)): # cluster must have at least 5% of the points within or be removed
                clusterIndex.append(i)
                redo = True


        # delete clusters if necessary
        if(redo):
            for i in range(len(clusterIndex) - 1, 0 - 1, -1):
                clusters.pop(clusterIndex[i])

            # reassign points
            for point in data:
                minValue = math.inf
                bestClus = -1
                for i in range(len(clusters)):
                    test = clusters[i].pointDistance(point)
                    if test < minValue:
                        minValue = test
                        bestClus = i
                point.assignCluster(bestClus)

    #Maximization Step
    def M_step(self, clusters, data):

        # Get new mu for each cluster
        for i in range(len(clusters)):
            counter = 0
            total = [0] * len(clusters[i].mu)
            for point in data:
                if(point.cluster == i):
                    counter = counter + 1
                    for j in range(len(point.features)):
                        total[j] = total[j] + point.features[j]
            for k in range(len(clusters[i].mu)):
                clusters[i].mu[k] = total[k] / counter

    #Trains data if there is a testSet present
    def Train(self, clusters, data, trainingSets):
        assignedClusters = []

        for set in trainingSets:
            set.cluster = -1
            minValue = math.inf
            bestClus = -1
            for i in range(len(clusters)):
                if i in assignedClusters:
                    continue
                test = set.clusterDistance(clusters[i])
                if test < minValue:
                    minValue = test
                    bestClus = i

            for index in set.points:
                data[index].cluster = bestClus

            set.cluster = bestClus
            assignedClusters.append(bestClus)

    #create testSet objects based on data
    def trainingSetsCreation(self, data, testSetIndexes, dataClassification):
        #find exisiting groups in classification
        groupContainer = []
        for group in dataClassification:
            if(group not in groupContainer):
                groupContainer.append(group)
        groupContainer.sort()


        #create testset objects
        testSets = []
        for group in groupContainer:
            count = 0
            mu = [0] * len(data[0])
            indexCollector = []
            for index in testSetIndexes:
                if(dataClassification[index] == group):
                    count = count + 1
                    indexCollector.append(index)
                    for i in range(len(data[index])):
                        mu[i] = mu[i] + data[index][i]

            if count != 0:
                for i in range(len(mu)):
                    mu[i] = mu[i] / count
                testSets.append(TestSet(group, indexCollector, mu))

        return testSets

    #create point objects based on data
    def pointCreation(self, data):
        points = []
        for row in data:
            points.append(Point(row))
        return points

    #create cluster objects based of data
    def clusterCreation(self, muValues): #create cluster objects based of data
        clusters = []
        for row in muValues:
            clusters.append(Cluster(row))
        return clusters

    #return objects for EM
    def info(self):
        totalClus = len(self.clusters)
        clusCenters = []
        for cluster in self.clusters:
            clusCenters.append(cluster.mu)
        classifications = []
        for point in self.dataPoints:
            classifications.append(point.cluster)
        return(totalClus, clusCenters, classifications)

    #convert to pandas dataframe
    def pandasConvert(self,data, clusterData):
        data = data.to_numpy().tolist()
        clusterData = clusterData.to_numpy().tolist()
        return(data, clusterData)
