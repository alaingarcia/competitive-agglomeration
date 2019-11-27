import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot(in_data, cluster_num1, center_list1, cluster_num2=None, center_list2=None, Dim=2, color=("blue","orange")):
    fig = plt.figure()
    ax = plt.axes()
    x = in_data.iloc[:,0]
    y = in_data.iloc[:,1]
    ax.scatter(x, y,s=1, color="black")
    
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