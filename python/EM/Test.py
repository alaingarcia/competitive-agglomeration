from datetime import datetime as dt
from sklearn import datasets #For iris dataset
import random
from Expectation_Maximization import ExpectationMaximization

if __name__ == '__main__':

    print("Creating testing data")
    data = datasets.load_iris().data #read dataset
    clusterCollect = []
    clusterCollect.append([4.9,2.8,5,2])
    clusterCollect.append([6,3.3,3.7,1])
    clusterCollect.append([4,2.9,4,1.])
    clusterCollect.append([4.2,4.1,5.7,1])
    clusterCollect.append([6,2,5,1])
    clusterCollect.append([6.5,3.2,4,1.1])

    #random test values
    testSetIndexes = []
    while(len(testSetIndexes) < (len(data)/3) ):
        randInt = random.randint(0, len(data) - 1)
        if(randInt not in testSetIndexes):
            testSetIndexes.append(randInt)

    testSetIndexes.sort()
    print(testSetIndexes)
    #end of random test values

    #random classification values
    testClassification = []
    for i in range(len(data)):
        testClassification.append(random.randint(0, 5))

    print(testClassification)

    print("End of testing data \n\n\n\n")

    startTime = dt.now()
    #EM = ExpectationMaximization(data,clusterCollect)
    EM = ExpectationMaximization(data,clusterCollect, testSet=testSetIndexes, classification=testClassification)
    endTime = dt.now()

    duration = endTime - startTime
    print("Time taken in milliseconds: {0}".format(duration.microseconds/1000))

    for cluster in EM.clusters:
        print(cluster.mu)

    for set in EM.trainingSets:
        print(set.cluster)
