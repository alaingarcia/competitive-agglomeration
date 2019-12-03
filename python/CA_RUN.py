from CA import CA
from accuracy import matrixAccuracy
import pandas as pd
import numpy as np
import time

def CA_RUN(cluster_num=100, iterations=50):
    # Read data (in_centers.csv, in_data.csv only things that are required)
    
    # A1 Synthetic Dataset Example
    #"""
    in_centers = pd.read_csv('data/A1-Dataset/random-centers-100.csv', header=None)
    in_data = pd.read_csv('data/A1-Dataset/a1.csv', header=None)
    actual_classification_file = pd.read_csv("data/A1-Dataset/a1-ga.csv", header=None)
    actual_classification_vector = actual_classification_file[0].tolist()
    actual_clusters = pd.read_csv("data/A1-Dataset/a1-ga-cb.csv", header=None).to_numpy().tolist()
    #"""
    
    # MATLAB Example (For this example, cluster num must be 10)
    """
    cluster_num = 10
    in_centers = pd.read_csv('data/MATLAB-Data-Example/InCenters.csv', header=None)
    in_data = pd.read_csv('data/MATLAB-Data-Example/InData.csv', header=None)
    actual_classification_file = pd.read_csv("data/MATLAB-Data-Example/ActualClassification.csv", header=None)
    actual_classification_vector = actual_classification_file[0].tolist()
    actual_clusters = pd.read_csv("data/MATLAB-Data-Example/OutCenters.csv", header=None)
    """
    
    # Calculate vector number, dimensions from shape
    vector_num = in_data.shape[0]
    dimensions = in_data.shape[1]
    in_centers = in_centers[0:cluster_num]

    # --------------------------CA-----------------------------------------
    print("Starting CA Algorithm:")
    CA_start_time = int(round(time.time() * 1000))
    CA_cluster_num, CA_center_list, CA_classification_list = CA.CA(in_data, in_centers, max_iterations=iterations,
                                                                    cluster_num=cluster_num, vector_num=vector_num,
                                                                    dimensions=dimensions)
    CA_end_time = int(round(time.time() * 1000))


    CA_time = (CA_end_time-CA_start_time)

    CA_accuracy = matrixAccuracy(actual_clusters, actual_classification_vector, CA_center_list, CA_classification_list)
    print("CA final cluster number: {}".format(CA_cluster_num))
    print("CA time: {} ms".format(CA_time))
    print("CA accuracy: {}%\n".format(CA_accuracy*100))

    """
    # Output to files
    pd.DataFrame(CA_center_list).to_csv('CA_center_list.csv', index=False)
    pd.DataFrame(CA_classification_list).to_csv('CA_classification_list.csv', index=False)

    all_data = pd.concat([in_data, actual_classification_file, pd.DataFrame(CA_classification_list)], axis=1)
    all_data.to_csv("all_data.csv")

    cluster_data = pd.concat([actual_clusters, pd.DataFrame(CA_center_list)], axis=1)
    cluster_data.to_csv("cluster_data.csv", index=False)
    """

    return(CA_cluster_num, CA_center_list)
    
if __name__ == '__main__':
    CA_RUN(cluster_num=30)