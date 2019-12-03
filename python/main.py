from CA_RUN import CA_RUN
from EM_RUN import EM_RUN
from plot import plot
import pandas as pd
import numpy as np
import time

if __name__ == '__main__':
    in_data = pd.read_csv('data/A1-Dataset/a1.csv', header=None)
    actual_clusters = pd.read_csv("data/A1-Dataset/a1-ga-cb.csv", header=None)
    # Graph the location of the actual clusters
    plot(in_data, actual_clusters=actual_clusters, save=True, save_name="result_plots/actual_plot.png", show=False)
    
    # Run EM and make a graph
    EM_cluster_num, EM_center_list = EM_RUN(NumClust=30)
    plot(in_data, cluster_num1=EM_cluster_num, center_list1=EM_center_list, save=True, save_name="result_plots/EM_plot.png", color=("orange",), show=False)
    
    # Run CA and make a graph
    CA_cluster_num, CA_center_list = CA_RUN(cluster_num=30)
    plot(in_data, cluster_num1=CA_cluster_num, center_list1=CA_center_list, save=True, save_name="result_plots/CA_plot.png", color=("red",), show=False)
   
    #"""
    # Plot everything together (actual cluster centers, EM Clusters, CA Clusters)
    plot(in_data, actual_clusters=actual_clusters, 
        cluster_num1=CA_cluster_num, center_list1=CA_center_list,
        cluster_num2=EM_cluster_num, center_list2=EM_center_list,
        save=True, save_name="result_plots/everything_plot.png", show=False)
    #"""

    print("Finished running both algorithms. Please check result_plots folder for plots.")

    """
    timing = dict()
    for i in range(1,100):
        timing[i] = []
        for j in range(0,3):
            timing[i].append(CA_RUN(cluster_num=i))

    df = pd.DataFrame.from_dict(timing)
    df.to_csv("timing.csv", index=False, header=None)
    """