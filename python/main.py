from CA_RUN import CA_RUN
from EM_RUN import EM_RUN
from plot import plot
import pandas as pd
import numpy as np
import time

if __name__ == '__main__':
    in_data = pd.read_csv('data/training/a1.csv', header=None)
    #CA_cluster_num, CA_center_list = CA_RUN()
    #EM_cluster_num, EM_center_list = EM_RUN()
    #plot(in_data, CA_cluster_num, CA_center_list, EM_cluster_num, EM_center_list)
    #plot(in_data, CA_cluster_num, CA_center_list)
    #plot(in_data)

    timing = dict()
    for i in range(1,100):
        timing[i] = []
        for j in range(0,3):
            timing[i].append(CA_RUN(cluster_num=i))

    df = pd.DataFrame.from_dict(timing)
    df.to_csv("timing.csv", index=False, header=None)