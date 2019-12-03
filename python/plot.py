import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot(in_data, actual_clusters=None, cluster_num1=None, center_list1=None, cluster_num2=None, center_list2=None, Dim=2, color=("red","orange"), save=False, save_name=None, show=True):
    fig = plt.figure()
    ax = plt.axes()
    x = in_data.iloc[:,0]
    y = in_data.iloc[:,1]
    ax.set_facecolor("#c4c4c4")
    ax.scatter(x, y,s=1, color="black")

    if actual_clusters is not None:
        x = actual_clusters.iloc[:,0]
        y = actual_clusters.iloc[:,1]
        ax.scatter(x, y, color="lime")
    
    if(cluster_num1 and center_list1):
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


    if show:
        """
        legend = ax.legend({"Lmao":"blue"} ,loc="upper left", title="Cluster Centers")
        legend.legendHandles[0]._legmarker.set_markersize(6)
        for label in legend.legendHandles:
            label._legmarker.set_markersize(6)
        ax.add_artist(legend)
        """
        plt.show()

    if save:
        if save_name is not None:
            plt.savefig(save_name)

if __name__ == "__main__":
    in_data = pd.read_csv("data/training/a1.csv", header=None)
    actual_clusters = pd.read_csv("data/training/a1-ga-cb.csv", header=None)
    plot(in_data, actual_clusters)