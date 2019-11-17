from CA import CA
from EM import EM
import pandas as pd


if __name__ == '__main__':
    # Read data (InCenters.csv, InData.csv only things that are required)
    InCenters = pd.read_csv('data/InCenters.csv', header=None)
    InData = pd.read_csv('data/InData.csv', header=None)
    #MinClSize = np.zeros((50,2))

    # Get results
    ExMa = EM.ExpectationMaximization(InData, InCenters, rep=50)
    NumClust, Center = CA.CA(InData, InCenters)
