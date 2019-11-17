from CA import CA
from EM import EM
import pandas as pd
import time


if __name__ == '__main__':
    # Read data (InCenters.csv, InData.csv only things that are required)
    InCenters = pd.read_csv('data/InCenters.csv', header=None)
    InData = pd.read_csv('data/InData.csv', header=None)

    # Get results and timing
    CA_start_time = time.time_ns()
    NumClust, OutCenter, Classifications = CA.CA(InData, InCenters, MaximumIt=50)
    CA_end_time = time.time_ns()

    # Print results
    print("CA final number of cluster: {}".format(NumClust))
    print("CA cluster centers: {}".format(OutCenter))
    print("CA classification vector: {}".format(Classifications))
    
    # Convert ns to ms
    CA_time = (CA_end_time-CA_start_time)/1e6

    #ExMa = EM.ExpectationMaximization(InData, InCenters, rep=50)

    def accuracy(actual_classification_vector, attempt_classification_vector):
        correct = 0
        for i in range(0, len(actual_classification_vector)):
            if actual_classification_vector[i] == attempt_classification_vector[i]:
                correct += 1
        return correct/len(actual_classification_vector)

    actual_classification_file = pd.read_csv("data/ActualClassification.csv", header=None)
    actual_classification_vector = actual_classification_file[0].tolist()
    CA_accuracy = accuracy(actual_classification_vector,Classifications)

    print("CA time: {} ms".format(CA_time))
    print("CA accuracy: {}%".format(CA_accuracy*100))