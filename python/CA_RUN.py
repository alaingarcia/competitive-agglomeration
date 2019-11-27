from CA import CA
from accuracy import matrixAccuracy
import pandas as pd
import numpy as np
import time

def CA_RUN(cluster_num=100, iterations=50):
    # Read data (in_centers.csv, in_data.csv only things that are required)
    in_centers = pd.read_csv('data/training/random-centers-100.csv', header=None)
    in_data = pd.read_csv('data/training/a1.csv', header=None)
    vector_num=in_data.shape[0]
    in_centers = in_centers[0:cluster_num]
    actual_classification_file = pd.read_csv("data/training/a1-ga.csv", header=None)
    actual_classification_vector = actual_classification_file[0].tolist()
    actual_clusters = pd.read_csv("data/training/a1-ga-cb.csv", header=None).to_numpy().tolist()

    # --------------------------CA-----------------------------------------
    print("Starting CA Algorithm:")
    CA_start_time = time.time_ns()
    CA_cluster_num, CA_center_list, CA_classification_list = CA.CA(in_data, in_centers, max_iterations=iterations, cluster_num=cluster_num, vector_num=vector_num)
    CA_end_time = time.time_ns()

    """
    print("CA final number of cluster: {}".format(CA_cluster_num))
    print("CA cluster centers: {}".format(CA_center_list))
    print("CA classification vector: {}".format(CA_classification_list))
    """
    # Convert ns to ms
    CA_time = (CA_end_time-CA_start_time)/1e6

    #CA_accuracy = matrixAccuracy(actual_clusters, actual_classification_vector, CA_center_list, CA_classification_list)
    print(cluster_num)
    print("CA time: {} ms".format(CA_time))
    #print("CA accuracy: {}%\n".format(CA_accuracy*100))

    """
    pd.DataFrame(CA_center_list).to_csv('CA_center_list.csv', index=False)
    pd.DataFrame(CA_classification_list).to_csv('CA_classification_list.csv', index=False)

    all_data = pd.concat([in_data, actual_classification_file, pd.DataFrame(CA_classification_list)], axis=1)
    all_data.to_csv("all_data.csv")

    cluster_data = pd.concat([actual_clusters, pd.DataFrame(CA_center_list)], axis=1)
    cluster_data.to_csv("cluster_data.csv", index=False)
    """
    #return(CA_cluster_num, CA_center_list)
    return(CA_time)
if __name__ == '__main__':
    CA_RUN()