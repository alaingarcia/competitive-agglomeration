from CA import CA
from EM import EM
import pandas as pd


if __name__ == '__main__':
    # Read data (InCenters.csv, InData.csv only things that are required)
    InCenters = pd.read_csv('data/InCenters.csv', header=None)
    InData = pd.read_csv('data/InData.csv', header=None)

    # Get results
    NumClust, OutCenter, Classifications = CA.CA(InData, InCenters, MaximumIt=50)
    print("CA final number of cluster: {}".format(NumClust))
    print("CA cluster centers: {}".format(OutCenter))
    print("CA classification vector: {}".format(Classifications))

    ExMa = EM.ExpectationMaximization(InData, InCenters, rep=50)