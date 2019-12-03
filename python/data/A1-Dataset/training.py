import pandas as pd
import math

cluster_centers = pd.read_csv("a1-ga-cb.csv", names=["x","y"])
data_points = pd.read_csv("a1.csv", names=["x","y"])
classifications = pd.read_csv("a1-ga.csv", names=["cluster"])

data_with_labels = pd.concat([data_points, classifications], axis=1)
smallest_points = pd.DataFrame()
index_series = pd.Series()

for index, row in cluster_centers.iterrows():
    distances = ((data_points["x"]-row["x"]).pow(2) + (data_points["y"]-row["y"]).pow(2)).pow(0.5)
    indeces = distances.nsmallest(15).index
    x = data_points["x"].iloc[indeces]
    y = data_points["y"].iloc[indeces]
    temp_smallest = pd.concat([x,y], axis=1)
    smallest_points = smallest_points.append(temp_smallest)
    index_series = index_series.append(pd.Series(indeces))

    # DEBUG PRINTS
    #print("Closest point to ({},{}) is".format(row["x"], row["y"]))
    #print(temp_smallest)

smallest_points.to_csv("training.csv", index=False, header=None)
index_series.to_csv("training_index.csv", index=False, header=None)