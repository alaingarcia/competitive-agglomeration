from CA import CA
from EM import Expectation_Maximization
import pandas as pd


if __name__ == '__main__':
    # Read data (InCenters.csv, InData.csv only things that are required)
    InCenters = pd.read_csv('data/InCenters.csv', header=None)
    InData = pd.read_csv('data/InData.csv', header=None)
    #MinClSize = np.zeros((50,2))

    # Get results
    NumClust, Center = CA.CA(InData, InCenters)
    print
    