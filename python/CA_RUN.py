from CA import CA
from accuracy import matrixAccuracy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time

def plot(in_data, cluster_num1, center_list1, cluster_num2=None, center_list2=None, Dim=2, color=("red","black")):
    fig = plt.figure()
    ax = plt.axes()
    x = in_data.iloc[:,0]
    y = in_data.iloc[:,1]
    ax.scatter(x, y,s=1)
    
    centerCoord1 = np.zeros(Dim)
    for i in range(0, cluster_num1):
        for j in range(0, Dim):
            centerCoord1[j] = center_list1[i][j]
        ax.scatter(centerCoord1[0], centerCoord1[1], color=color[0])
    
    if(cluster_num2 and center_list2):
        centerCoord2 = np.zeros(Dim)
        for i in range(0, cluster_num2):
            for j in range(0, Dim):
                centerCoord2[j] = center_list2[i][j]
            ax.scatter(centerCoord2[0], centerCoord2[1], color=color[1])
    plt.show()

if __name__ == '__main__':
    # Read data (in_centers.csv, in_data.csv only things that are required)
    in_centers = pd.read_csv('data/training/random-centers-40.csv', header=None)
    in_data = pd.read_csv('data/training/a1.csv', header=None)
    cluster_num=25
    iterations=50
    vector_num=3000
    in_centers = in_centers[0:cluster_num]
    actual_classification_file = pd.read_csv("data/training/a1-ga.csv", header=None)
    actual_classification_vector = actual_classification_file[0].tolist()
    actual_clusters = pd.read_csv("data/training/a1-ga-cb.csv", header=None)

    # --------------------------CA-----------------------------------------
    print("Starting CA Algorithm:")
    CA_start_time = time.time_ns()
    CA_cluster_num, CA_center_list, CA_classification_list = CA.CA(in_data, in_centers, max_iterations=iterations, cluster_num=cluster_num, vector_num=vector_num)
    CA_end_time = time.time_ns()

    # Print results
    print("CA final number of cluster: {}".format(CA_cluster_num))
    print("CA cluster centers: {}".format(CA_center_list))
    print("CA classification vector: {}".format(CA_classification_list))

    # Convert ns to ms
    CA_time = (CA_end_time-CA_start_time)/1e6

    #pd.DataFrame(CA_center_list).to_csv('CA_center_list.csv', index=False)
    #pd.DataFrame(CA_classification_list).to_csv('CA_classification_list.csv', index=False)

    #CA_accuracy = matrixAccuracy(actual_clusters, actual_classV, attempt_clusters, attempt_classV)
    print("CA time: {} ms".format(CA_time))
    #print("CA accuracy: {}%\n".format(CA_accuracy*100))

    all_data = pd.concat([in_data, actual_classification_file, pd.DataFrame(CA_classification_list)], axis=1)
    all_data.to_csv("all_data.csv")

    cluster_data = pd.concat([actual_clusters, pd.DataFrame(CA_center_list)], axis=1)
    cluster_data.to_csv("cluster_data.csv", index=False)

    plot(in_data, CA_cluster_num, CA_center_list)