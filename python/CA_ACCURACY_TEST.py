import pandas as pd
from CA import CA
from accuracy import matrixAccuracy

all_data = pd.read_csv("all_data.csv")
cluster_data = pd.read_csv("cluster_data.csv")

actual_clusters = list()
attempt_clusters = list()
for index, row in cluster_data.iterrows():
    if not pd.isna(row["x_true"]) and not pd.isna(row["y_true"]):
        actual_clusters.append([row["x_true"],row["y_true"]])
    if not pd.isna(row["x_pred"]) and not pd.isna(row["y_pred"]):
        attempt_clusters.append([row["x_pred"],row["y_pred"]])

attempt_classV = all_data["pred"].tolist()
actual_classV = all_data["TRUE"].tolist()

CA_accuracy = matrixAccuracy(actual_clusters, actual_classV, attempt_clusters, attempt_classV)
print(CA_accuracy)