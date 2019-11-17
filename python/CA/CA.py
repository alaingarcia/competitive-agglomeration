import numpy as np
import math
import matplotlib.pyplot as plt

EPS = 1.0e-6
MAX_CENT_DIFF = 0.001
INFINITY = 1.0e35

class Feature_Info():
    def __init__(self, memship, dist, dimen, cluster):
        self.memship = memship
        self.dist = dist
        self.dimen = dimen
        self.cluster = cluster

"""
/**************************************************************************************************/
/***  This file contains procedures that perform fuzzy clustering with aglomeration (CA).		***/
/***  These include:  CA.							***/
/***         																					***/
/***  		Author:               Hichem FRIGUI													***/
/***        Date last updated:    September 28th 2005											***/
/**************************************************************************************************/


/**************************************************************************************************/
/***  This procedure implements the CA algorithm.												***/
/***																							***/
/***	INPUT:																					***/
/***	   FeatVect:  A structure that holds all feature vectors. 								***/
/***	   NumClust:  Total number of clusters.													***/
/***	   m:	      fuzzifier.																***/
/***	   Center:    Initial set of centers.  If not included centers will be initialized.		***/
/***	   MaxIter:   Maximum number of iteration.												***/
/***																							***/
/***	OUTPUT:																					***/
/***	   Center:   Final centers of the clusters.												***/
/***	   FeatVect: This structure contains the final partition and							***/
/***		     the membership assignment.														***/
/**************************************************************************************************/
"""
def CA(InData, InCenters,
        MaximumIt=50, Eita_0=2.5, TAU=10, EXPINC=25, 
        EXPDEC=35, EPS=1.0e-6, MAX_CENT_DIFF=0.001,
        NumClust=10, NumOfVectors=60, Dim=2, m_p=2):
    
    # Initialize Data (from InData)
    FeatVect = []
    for i in range(0, NumOfVectors):
        FeatVect.insert(i,Feature_Info(np.zeros(NumClust), 
                        np.zeros(NumClust), np.zeros(Dim), None))
    # Populate
    for i in range(0, Dim):
        for j in range(0, NumOfVectors):
            FeatVect[j].dimen[i] = InData.iloc[j,i]

    # Initialize Centers (random data from InCenters)
    Center = InCenters.values

    # MinPts are zero vectors ?
    MinPts = np.zeros(MaximumIt)
    MinPts2 = np.ones(MaximumIt)
    MinPts2[int(MaximumIt/2):MaximumIt] *= 5

    i, j, IterNum = 0, 0, 0
    AggCte = 0.0
    PreviousNumClust = NumClust

	# PreviousCenter = DMatrix_2D(0,NumClust,0,Dim)
    PreviousCenter = np.zeros((NumClust, Dim))

    for i in range(0, NumClust):
        for j in range(0, Dim):
            PreviousCenter[i][j] = Center[i][j]

    while (IterNum < MaximumIt):
        #print("IterNum={}\n".format(IterNum))
        
        FeatVector, Center = EucDistance(FeatVect, Center, NumClust, NumOfVectors, Dim, m_p)
        FeatVector = AssignPts(FeatVect, NumClust, NumOfVectors)

        if (IterNum>2):
            FeatVector, AggCte = Get_AggCte(FeatVect,NumClust,IterNum,NumOfVectors,Eita_0, TAU, EXPINC, EXPDEC)
            #print("AggConst={}\n".format(AggCte))
            FeatVector = CompMem(FeatVect, NumClust,AggCte,NumOfVectors)
            FeatVector, NumClust, MinPts, MinPts2 = UpdateNumClusters(FeatVect, NumClust, IterNum, NumOfVectors, MinPts, MinPts2)
        else:
            FeatVector = FuzzMem(FeatVect, NumClust, NumOfVectors)
        
        FeatVector, Center = FuzzCenters(FeatVect,Center,NumClust,NumOfVectors,Dim)
        
        CenDiff = 0.0
        for i in range(0, NumClust):
            for j in range(0, Dim):
                CenDiff += abs(PreviousCenter[i][j]-Center[i][j])
                
        if (CenDiff < MAX_CENT_DIFF and PreviousNumClust == NumClust):
            break
        
        IterNum += 1 
        #print("[It# {}]Num. Clust={}\n".format(IterNum,NumClust))
        PreviousNumClust = NumClust

		# Copy the centers into the previous centers
        for i in range(0, NumClust):
            for j in range(0, Dim):
                PreviousCenter[i][j] = Center[i][j]
    
    FeatVector, Center = EucDistance(FeatVect, Center, NumClust, NumOfVectors, Dim, m_p)
    FeatVector = AssignPts(FeatVect, NumClust, NumOfVectors)
    
    Classifications = []
    for i in range(0, NumOfVectors):
        Classifications.append(FeatVect[i].cluster)

    OutCenters = []
    for i in range(0, NumClust):
        OutCenters.append([])
        for j in range(0, Dim):
            OutCenters[i].append(Center[i][j])
    """
    # ---------- OUTPUT TO FILES ---------------
    print(NumClust, file=open('NumClust.txt', 'w'))

    clusterOut = open('OutCluster.txt', 'w')
    for i in range(0, NumOfVectors):
        clusterOut.write(str(FeatVect[i].cluster) + '\n')
    clusterOut.close()
    
    centerOut = open('OutCenters.txt', 'w')
    for i in range(0, NumClust):
        for j in range(0, Dim):
            centerOut.write(str(Center[i][j]) + ' ')
        centerOut.write('\n')
    centerOut.close()

    # ------------- PLOT --------------
    
    fig = plt.figure()
    ax = plt.axes()
    x = InData.iloc[:,0]
    y = InData.iloc[:,1]
    ax.scatter(x, y)
    
    centerCoord = np.zeros(Dim)
    for i in range(0, NumClust):
        for j in range(0, Dim):
            centerCoord[j] = Center[i][j]
        ax.scatter(centerCoord[0], centerCoord[1], color='red')
    plt.show()
    """

    return(NumClust, OutCenters, Classifications)
"""
/**************************************************************************************************/
/*** This procedure computes the Euclidean distance of all feature vectors						***/
/*** from all clusters. Then, it crisply assigns each feature point to the						***/
/*** closest cluster.																			***/
/***																							***/
/***     INPUT:																					***/
/***         FeatVector:    A structure containing the feature vectors.							***/
/***         NumClust:      Total number of clusters.											***/
/***         Center:        Centers of all clusters.  											***/
/***		 NumOfVectors:  Total number of feature vectors.									***/
/***		 Dim:			Total number of dimensions.											***/
/***																							***/
/***     OUTPUT:																				***/
/***         FeatVector:    The component of this structure that holds the						***/
/***                        distance, and the cluster assignment.								***/
/**************************************************************************************************/
"""

def EucDistance(FeatVector, Center, NumClust, NumOfVectors, Dim, m_p):
    for i in range(0, NumOfVectors):
        for j in range(0, NumClust):
            FeatVector[i].dist[j] = 0
            for k in range(0, Dim):
                temp = FeatVector[i].dimen[k]-Center[j][k]
                temp1 = 1
                for p in range(0, m_p):
                    temp1 = temp1*temp
                FeatVector[i].dist[j] += abs(temp1)
            FeatVector[i].dist[j] = max(EPS, FeatVector[i].dist[j])
    return(FeatVector, Center)

"""
/**************************************************************************************************/
/*** This procedure crisply assigns each feature point to the closest cluster					***/
/***																							***/
/***     INPUT:																					***/
/***         FeatVector:    A structure containing the feature vectors.							***/
/***         NumOfVectors:  Total number of feature vectors.									***/
/***         NumClust: Total number of clusters.												***/
/***     OUTPUT:																				***/
/***         FeatVector:    The component of this structure that holds the						***/
/***                        cluster assignment.             									***/
/**************************************************************************************************/
"""

def AssignPts(FeatVector, NumClust, NumOfVectors):
    """
    int i, j, index; double min;
	for (i=0; i<NumOfVectors; i++){
		for (j=0, min=INFINITY; j<NumClust; j++){
			if ( (FeatVector+i)->dist[j] <  min ) {
				min = (FeatVector+i)->dist[j]; 
				index = j; 
			}
			(FeatVector+i)->cluster = index;
		}
	}
    """
    for i in range(0, NumOfVectors):
        min = INFINITY
        for j in range(0, NumClust):
            if FeatVector[i].dist[j] < min:
                min = FeatVector[i].dist[j]
                index = j
            FeatVector[i].cluster = index
    return(FeatVector)

"""
/**************************************************************************************************/
/*** This function computes the agglomeration constant. This constant is used					***/
/*** to control the cometition rate in the competitive clustering Alg.							***/
/***									    													***/
/***      INPUT:								  												***/
/*** 	     FeatVector:   A structure that has all feature vectors.        					***/
/***         NumOfVectors: Total number of feature vectors.		    							***/
/***         NumOfClusters:Total number of clusters.			    							***/
/***         IterNum:      Current iteration number.			    							***/
/***         Card:         Cardinalities of all clusters.                    					***/ 
/***         flg:          0 for non robust and 1 for robust clustering      					***/
/*** 									   														***/
/***      RETURNS:							    												***/
/***         The agglomeration constant.                                     					***/
/**************************************************************************************************/
"""
def Get_AggCte (FeatVector, NumClust, IterNum, NumOfVectors, Eita_0, TAU, EXPINC, EXPDEC):
    ObjFunc, SumCard2 = 0, 0
    Card = np.zeros(NumClust)

    for i in range(0, NumClust):
        Card[i] = 0
        for j in range(0, NumOfVectors):
            temp = FeatVector[j].memship[i]
            Card[i] += temp
            ObjFunc += temp * temp * FeatVector[j].dist[i]
        SumCard2 += Card[i] * Card[i]

    if (IterNum < EXPINC):
        Exponent = (IterNum - EXPINC)/TAU  
    elif (IterNum < EXPDEC):
        Exponent = 0
    else:
        Exponent = (EXPDEC - IterNum) / TAU
    
    alpha = Eita_0 * math.exp(Exponent) * ObjFunc/SumCard2
    
    return(FeatVector, alpha)

"""
/**************************************************************************************************/
/*** This procedure computes the membership of all feature vectors in all						***/
/*** clusters. The computed mebership has 2 terms: the Fuzzy membership and   					***/
/*** a bias term.                                                             					***/
/***                                                                          					***/
/***     INPUT:								    												***/
/***         FeatVector:    A structure containing the feature vectors.       					***/
/***         NumOfVectors:  Total number of feature vectors.                  					***/
/***         NumOfClusters: Total number of clusters.			    							***/
/***         AggCte:        A constant term that controls the rate of         					***/
/***                        agglomeration.                                    					***/
/***     OUTPUT:								    											***/
/***         FeatVector:    The membership components of this structure are   					***/
/***                        updated.                                          					***/
/**************************************************************************************************/
"""
def CompMem(FeatVector, NumClust, AggCte, NumOfVectors):
    Card = np.zeros(NumClust)
    
    for j in range(0, NumClust):
        Card[j] = 0
        if (AggCte >= 0.0):
            for i in range(0, NumOfVectors):
                Card[j] += FeatVector[i].memship[j]
    
    for i in range(0, NumOfVectors):
        Max, Min, SumInvDist, AvgCard = 1.0, 0.0, 0.0, 0.0
        for j in range(0, NumClust):
            SumInvDist += 1.0/FeatVector[i].dist[j]
        for j in range(0, NumClust):
            AvgCard += Card[j]/FeatVector[i].dist[j]
        AvgCard /= SumInvDist

        for j in range(0, NumClust):
            Mem_FCM = 1.0 / (SumInvDist * FeatVector[i].dist[j])
            Mem_Bias = (Card[j] - AvgCard) / FeatVector[i].dist[j]
            FeatVector[i].memship[j] = Mem_FCM + AggCte * Mem_Bias

            # Force membership to be in the interval [0, 1]
            FeatVector[i].memship[j] = max(0, FeatVector[i].memship[j])
            FeatVector[i].memship[j] = min(1, FeatVector[i].memship[j])

        SumMem = 0
        for j in range(0, NumClust):
            SumMem += FeatVector[i].memship[j]

        for j in range(0, NumClust):
            FeatVector[i].memship /= SumMem
        
        SumMem = 0
        for j in range(0, NumClust):
            SumMem += FeatVector[i].memship[j]
        
        if SumMem > 1.0e-10:
            for j in range(0, NumClust):
                FeatVector[i].memship[j] /= SumMem
        else:
            for j in range(0, NumClust):
                FeatVector[i].memship[j] = 1 / NumClust
    return(FeatVector)
"""
/**************************************************************************************************/
/*** This procedure updates the number of clusters for the competitive alg.						***/
/*** A cluster is deleted if its cardinality is less than its dimensionality. 					***/
/***                                                                          					***/
/***     INPUT:                                                               					***/
/***       FeatVect:    A structure holding the feature vectors.              					***/
/***       NumOfVectors:     Total number of feature vectors.                      				***/
/***       NumClust:    Initial number of clusters.                           					***/
/***       Dim:         Dimensionality of the feature space.                  					***/
/***                                                                          					***/
/***     OUTPUT:                                                              					***/
/***       NumClust:    Updated number of clusters.                           					***/
/**************************************************************************************************/
"""
def UpdateNumClusters(FeatVect, NumClust, IterNum, NumOfVectors, MinPts, MinPts2):
    i, j, OptNumClust = (0, 0, 0)
    Empty_clust = np.zeros(NumClust)
    Card = np.zeros(NumClust)
    Card2 = np.zeros(NumClust)
    Card3 = np.zeros(NumClust)
    
    
    for i in range(0, NumClust):
        Empty_clust[i] = 0 # Assume first that all clusters are not empty
        Card[i] = 0
        Card3[i] = 0
        for j in range(0, NumOfVectors):
            Card[i] += FeatVect[j].memship[i]*FeatVect[j].memship[i]
            Card3[i] += FeatVect[j].memship[i]
            if FeatVect[j].cluster == i:
                Card2[i] += 1
        
        if (Card[i] <= MinPts[IterNum]) or (Card2[i] <= MinPts2[IterNum]):
            Empty_clust[i] = 1

        #if IterNum >= 3:
            #print('{} ===> Card[{}] = {}   Card2[{}] = {}\n'.format(IterNum, i, Card[i], i, Card2[i]))

    for i in range(0, NumClust):
        if Empty_clust[i] == 0:
            for j in range(0, NumOfVectors):
                FeatVect[j].memship[OptNumClust] = FeatVect[j].memship[i]
                if FeatVect[j].cluster == i:
                    FeatVect[j].cluster = OptNumClust
            OptNumClust += 1

    return(FeatVect, OptNumClust, MinPts, MinPts2)

"""
/**************************************************************************************************/
/*** This procedure computes the fuzzy membership of all feature vectors in						***/
/*** all clusters.  						           											***/
/***                                                                          					***/
/***     INPUT:								    												***/
/***         FeatVector:    A structure containing the feature vectors.       					***/
/***        NumOfVectors:  Total number of feature vectors.                  					***/
/***         NumOfClusters: Total number of clusters.			    							***/
/***         m:             The fuzzifier.				    									***/
/***     OUTPUT:								    											***/
/***         FeatVector:    The membership components of this structure are   					***/
/***                        updated.                                          					***/
/**************************************************************************************************/
"""

def FuzzMem(FeatVector, NumClust, NumOfVectors):
    for i in range(0, NumOfVectors):
        SumInvDist = 0
        for j in range(0, NumClust):
            SumInvDist +=  1 / FeatVector[i].dist[j]
        for j in range(0, NumClust):
            temp = FeatVector[i].dist[j]
            FeatVector[i].memship[j] = 1 / (SumInvDist*temp)
    return(FeatVector)
"""
/**************************************************************************************************/
/*** This procedure updates the centers for the Fuzzy algorithms.								***/
/***                                                                          					***/
/***     INPUT:                                                               					***/
/***       FeatVect:    A structure holding the feature vectors.              					***/
/***       NumOfVectors:     Total number of feature vectors.                    				***/
/***       NumOfClusters:    Total number of clusters.                             				***/
/***       Dim:         Dimensionality of the feature space.                  					***/
/***       m:           Fuzzifier.			                    								***/
/***                                                                          					***/
/***     OUTPUT:                                                              					***/
/***      Center: Updated centers.                          		    						***/
/**************************************************************************************************/
"""

def FuzzCenters(FeatVect, Center, NumClust, NumOfVectors, Dim):
    Card = np.zeros(NumClust)

    for j in range(0, NumClust):
        Card[j] = 0
        for k in range(0, Dim):
            Center[j][k] = 0

        for i in range(0, NumOfVectors):
            temp = FeatVect[i].memship[j] * FeatVect[i].memship[j]

            Card[j] += temp
            for k in range(0, Dim):
                Center[j][k] += temp * FeatVect[i].dimen[k]

        for k in range(0, Dim):
            Center[j][k] /= Card[j] + EPS

    return (FeatVect, Center)

"""
/**************************************************************************************************/
/*** This procedure initializes the centers for the C clusters. The centers						***/
/*** are chosen as C data points equally spaced.                              					***/
/***                                                                          					***/
/***     INPUT:                                                               					***/
/***       FeatVect:    A structure holding the feature vectors.              					***/
/***       NumOfVectors:     Total number of feature vectors.                 					***/
/***       NumOfClusters:    Total number of clusters.                        					***/
/***       Dim:         Dimensionality of the feature space.                  					***/
/***                                                                          					***/
/***     OUTPUT:                                                              					***/
/***       Center: Initial centers for the clusters.                          					***/
/**************************************************************************************************/
"""
def InitCenters (Center, FeatVect, NumClust, NumOfVectors, Dim):
    for k in range(0,Dim):
        for j in range(0, NumClust):
            index = j * NumOfVectors/NumClust
            Center[j][k] = FeatVect[index].dimen[k]
    return(Center, FeatVect)