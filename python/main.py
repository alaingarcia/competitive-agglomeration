from CA_RUN import CA_RUN
from EM_RUN import EM_RUN
from plot import plot
import pandas as pd
import numpy as np
import time

if __name__ == '__main__':
    in_data = pd.read_csv('data/training/a1.csv', header=None)
    CA_cluster_num, CA_center_list = CA_RUN()
    EM_cluster_num, EM_center_list = EM_RUN()
    plot(in_data, CA_cluster_num, CA_center_list, EM_cluster_num, EM_center_list)

