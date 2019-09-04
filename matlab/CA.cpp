#include <stdio.h>
#include <stdlib.h>

/**************************************************************************************************/
/***  This file contains procedures that perform fuzzy clustering with aglomeration (CA).		***/
/***  These include:  CA.							***/
/***         																					***/
/***  		Author:               Hichem FRIGUI													***/
/***        Date last updated:    September 28th 2005											***/
/**************************************************************************************************/

#include "CA.h"

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
int CA (struct Feature_Info *FeatVect, int NumOfVectors, int Dim, int NumClust , double **Center, int MaximumIt, double Eita_0, double *MinPts, double *MinPts2, double TAU, int EXPINC, int EXPDEC, double m_p)
{
	int i, j, IterNum = 0; 
	int PreviousNumClust;
	double AggCte,*Card,CenDiff;
	double **PreviousCenter;
	void InitCenters (double **,struct Feature_Info *,int,int,int);
	void EucDistance (struct Feature_Info *,double **,int,int,int,double);
	void AssignPts(struct Feature_Info *,int,int);
	double Get_AggCte (struct Feature_Info *,int,int,int,double, double, int, int);
	void CompMem (struct  Feature_Info *,int,double,int);
	void UpdateNumClusters (struct Feature_Info *,int *,int,int,double *,double *);
	void FuzzMem (struct Feature_Info *,int,int);
	void FuzzCenters (struct Feature_Info *,double **,int,int,int);
	double abs(double);

	AggCte=0.0;
	PreviousNumClust = NumClust;

	Card=DMatrix_1D(0,NumClust);

	PreviousCenter = DMatrix_2D(0,NumClust,0,Dim);
	/** Copy the centers into the previous centers**/
	for(i=0;i<NumClust;i++)
		for(j=0;j<Dim;j++) PreviousCenter[i][j] = Center[i][j];

	// Write output Centers
	//FILE *fp=fopen("OutCentersInit.txt","w");
	//for (j=0;j<NumClust;j++){
	//	for (i=0;i<Dim;i++){
 //			fprintf(fp,"%f ",Center[j][i]);
	//	}
	//	fprintf(fp,"\n");
	//}
	//fclose(fp);
    while (IterNum < MaximumIt){
        printf("IterNum=%d\n",IterNum);

		EucDistance (FeatVect,Center,NumClust,NumOfVectors,Dim, m_p);
		AssignPts   (FeatVect,NumClust,NumOfVectors);
      
		if (IterNum>2){
			AggCte = Get_AggCte (FeatVect,NumClust,IterNum,NumOfVectors,Eita_0, TAU, EXPINC, EXPDEC);
			printf("AggConst=%f\n",AggCte);
			CompMem (FeatVect, NumClust,AggCte,NumOfVectors);
			UpdateNumClusters (FeatVect,&NumClust,IterNum,NumOfVectors, MinPts, MinPts2);
		}
		else
			FuzzMem(FeatVect,NumClust,NumOfVectors);
   
		FuzzCenters (FeatVect,Center,NumClust,NumOfVectors,Dim);
   
		/** Check for convergence **/
		for (i=0, CenDiff=0.0; i<NumClust; i++)
			for (j=0; j<Dim; j++) CenDiff += abs(PreviousCenter[i][j]-Center[i][j]);
		if (CenDiff < MAX_CENT_DIFF && PreviousNumClust == NumClust) break;

		IterNum ++; 
		printf("[It# %d]Num. Clust=%d\n",IterNum,NumClust);
		PreviousNumClust = NumClust;
		/** Copy the centers into the previous centers**/
		for(i=0;i<NumClust;i++)
			for(j=0;j<Dim;j++) PreviousCenter[i][j] = Center[i][j];
		//if(IterNum==40)
		//{
		//	fp=fopen("OutCenters5.txt","w");
		//	for (j=0;j<NumClust;j++){
		//		for (i=0;i<Dim;i++){
 	//				fprintf(fp,"%f ",Center[j][i]);
		//		}
		//		fprintf(fp,"\n");
		//	}
		//	fclose(fp);
		//}
		if(IterNum==40){
		
		}	
	}
	EucDistance (FeatVect,Center,NumClust,NumOfVectors,Dim, m_p);
	AssignPts   (FeatVect,NumClust,NumOfVectors);
	
	return(NumClust);
}

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
void EucDistance (struct Feature_Info *FeatVector, double **Center, int NumClust, int NumOfVectors, int Dim, double m_p)
{
	int i, j, k; 
	double temp, temp1;

	for (i=0; i<NumOfVectors; i++){
		for (j=0; j<NumClust; j++){ 
			for (k=0,(FeatVector+i)->dist[j]=0.0; k<Dim; k++){
				temp = (FeatVector+i)->dimen[k]-Center[j][k];
				temp1 = 1;
				for (int p=0; p<m_p; p++){ 
				   temp1 = temp1*temp;
				}
				(FeatVector+i)->dist[j] += fabs(temp1);	  
			}
			(FeatVector+i)->dist[j] = MAX (EPS, (FeatVector+i)->dist[j]);
		}
	}
}

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
void AssignPts (struct Feature_Info *FeatVector, int NumClust, int NumOfVectors)
{
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
}

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
double Get_AggCte (struct Feature_Info *FeatVector, int NumClust, int IterNum, int NumOfVectors, double Eita_0, double TAU, int EXPINC, int EXPDEC)
{
    int i, j;
	double ObjFunc=0.0, SumCard2=0.0, alpha, Exponent; 
	double *Card = DMatrix_1D(0, NumClust);
	double exp(double);

	for (i=0; i<NumClust; i++){
		for (j=0, Card[i]=0; j<NumOfVectors; j++){
 			Card[i] += FeatVector[j].memship[i];
			ObjFunc += FeatVector[j].memship[i]*FeatVector[j].memship[i]*FeatVector[j].dist[i];
		}
		SumCard2 += Card[i]*Card[i];
	}

	if (IterNum < EXPINC){    
		Exponent = (IterNum-EXPINC)/TAU;
	}
	else if (IterNum <EXPDEC){
		Exponent = 0.0;
	}
	else{
		Exponent = (EXPDEC-IterNum)/TAU;
	}
	alpha = Eita_0  * exp(Exponent)  * ObjFunc/SumCard2;
	Free_DMatrix_1D(Card,0,NumClust);
	//return (alpha);
	return (alpha);
}

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

void CompMem (struct  Feature_Info *FeatVector, int NumClust, double AggCte, int NumOfVectors)
{
    int    i, j;
	double SumInvDist, AvgCard, Mem_FCM, Mem_Bias, SumMem;  
	double TotCard=0.0;
	double Max, Min;
	double *Card = DMatrix_1D(0, NumClust);

	/*** Compute Cardinality of each cluster ***/
	for (j=0; j<NumClust; j++){
		Card[j] = 0.0;
		if (AggCte>=0.0) 
			for (i=0; i<NumOfVectors; i++)      
				Card[j] +=  (FeatVector+i)->memship[j];
	}

  
	for (i=0; i<NumOfVectors; i++){
		Max = 1.0;  Min = 0.0;

		SumInvDist=0.0;   AvgCard=0.0;
		for (j=0; j<NumClust; j++)  
			SumInvDist += 1.0/(FeatVector+i)->dist[j];
		for (j=0; j<NumClust; j++)  
			AvgCard += Card[j]/(FeatVector+i)->dist[j];
		AvgCard /= SumInvDist;

		for (j=0; j<NumClust; j++) {
			Mem_FCM  = 1.0 / (SumInvDist*(FeatVector+i)->dist[j]);
			Mem_Bias = (Card[j]-AvgCard) / (FeatVector+i)->dist[j];
			(FeatVector+i)->memship[j] = Mem_FCM + AggCte*Mem_Bias; 

			/*** Force membership to be in the interval [0, 1] ***/
					//opt 1	
			         (FeatVector+i)->memship[j] = MAX (0.0, (FeatVector+i)->memship[j]);
			         (FeatVector+i)->memship[j] = MIN (1.0, (FeatVector+i)->memship[j]);
			
			//opt 2
			//if ( (FeatVector+i)->memship[j]>Max)    
			//	Max = (FeatVector+i)->memship[j];
			//if ( (FeatVector+i)->memship[j]<Min)    
			//	Min = (FeatVector+i)->memship[j];
		}
    
		for (j=0, SumMem=0.0; j<NumClust; j++) {
			//opt 2
			//(FeatVector+i)->memship[j] = ((FeatVector+i)->memship[j] - Min)/(Max-Min); 
			SumMem += (FeatVector+i)->memship[j];
		}
		for (j=0; j<NumClust; j++)   
			(FeatVector+i)->memship[j] /= SumMem;    

		for (j=0, SumMem=0.0; j<NumClust; j++)  
			SumMem += (FeatVector+i)->memship[j];
		if (SumMem > 1.0e-10)
			for (j=0; j<NumClust; j++)  
				(FeatVector+i)->memship[j] /= SumMem;
		else 
			for (j=0; j<NumClust; j++)  
				(FeatVector+i)->memship[j] = 1.0/NumClust; 
	}
	Free_DMatrix_1D(Card,0,NumClust);
}

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
void UpdateNumClusters (struct Feature_Info *FeatVect, int * NumClust, int IterNum, int NumOfVectors, double *MinPts, double *MinPts2)
{
	int i, j, OptNumClust=0, k;
	int *Empty_clust	= IMatrix_1D (0, (*NumClust)-1 );
	double *Card		= DMatrix_1D (0, (*NumClust)-1 ); 	
	int *Card2			= IMatrix_1D (0, (*NumClust)-1 ); 

	double *Card3		= DMatrix_1D (0, (*NumClust)-1 ); 

 	for (j=0;j<(*NumClust);j++)
		Card2[j]=0;

	for (i=0; i<(*NumClust); i++) {
		Empty_clust[i] = 0;   /** Assume first that all clusters are not empty **/
		for (j=0, Card[i]=0, Card3[i]=0 ; j<NumOfVectors; j++) {
			Card[i] += (FeatVect+j)->memship[i]*(FeatVect+j)->memship[i];
			Card3[i]+= (FeatVect+j)->memship[i]; 
			if ((FeatVect+j)->cluster==i)
				Card2[i]++;
		}
		if (Card[i]<MinPts[IterNum] || Card2[i]<MinPts2[IterNum])
			Empty_clust[i] = 1;

		if (IterNum>=3)
			printf ("%d ===> Card[%d] = %f   Card2[%d] = %d\n", IterNum, i, Card[i], i, Card2[i]); 
	}

//	//cheking for overlapping clusters
////	double SimClTh = 0.03;
//	//if ((IterNum>=25)&&((IterNum%5) == 0)){
//		for (i=0; i<(*NumClust); i++) {
//			for (j=i+1; j<(*NumClust); j++) {	
//				if ((Empty_clust[i] == 0)&&(Empty_clust[j] == 0)){
//					double t = 0;
//					for (k=0; k<NumOfVectors; k++) {
//						t += fabs((FeatVect+k)->memship[i] - (FeatVect+k)->memship[j]);
//					}			
//					t = t / (Card3[i] + Card3[j]);
//					if (t < SimClTh){
//						Empty_clust[j] = 1;
//						printf ("Merged ===> %d and  %d\n", i, j); 
//					}				
//				}
//			}
//		}
//	//}

	for (i=0; i<(*NumClust); i++) {
		/** Save parameters of non-empty clusters in consecutive order **/
		if (Empty_clust[i] == 0) {
			for (j=0; j<NumOfVectors; j++) {
				(FeatVect+j)->memship[OptNumClust] = (FeatVect+j)->memship[i];
				if ( (FeatVect+j)->cluster == i)    
					(FeatVect+j)->cluster = OptNumClust;
			}
			OptNumClust ++;
		}
	}
	Free_IMatrix_1D (Empty_clust, 0, (*NumClust)-1);
	Free_DMatrix_1D (Card, 0, (*NumClust)-1);
	Free_IMatrix_1D (Card2, 0, (*NumClust)-1);

	(*NumClust) = OptNumClust;	
}

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
void FuzzMem (struct Feature_Info *FeatVector, int NumClust, int NumOfVectors)
{
	int i, j;  
	double SumInvDist, temp;

	for (i=0; i<NumOfVectors; i++){
		for (j=0, SumInvDist=0.0; j<NumClust; j++)
			SumInvDist += 1.0/(FeatVector+i)->dist[j];
		for (j=0; j<NumClust; j++){
			temp = (FeatVector+i)->dist[j];
			(FeatVector+i)->memship[j] = 1.0/(SumInvDist*temp); 
		}
	}
}

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
void FuzzCenters (struct Feature_Info *FeatVect, double ** Center, int NumClust, int NumOfVectors, int Dim)
{
	int    i, j, k; 
	double temp;
	double *Card = DMatrix_1D(0,NumClust);

 
	for (j=0; j<NumClust; j++){
		Card[j] = 0.0;
		for (k=0; k<Dim; k++) 
			Center[j][k]=0.0;
 
		for (i=0; i<NumOfVectors; i++) {
			temp = (FeatVect+i)->memship[j]*(FeatVect+i)->memship[j];

			Card[j] += temp;
 			for (k=0; k<Dim; k++) 
				Center[j][k] += temp*(FeatVect+i)->dimen[k];
		}  
      
		for (k=0; k<Dim; k++)
			Center[j][k] /= (Card[j]+EPS);
	}
	Free_DMatrix_1D(Card,0,NumClust);
}

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
void InitCenters (double **Center, struct Feature_Info *FeatVect, int NumClust, int NumOfVectors, int Dim)
{
	int j, k,index;
    for (k=0;k<Dim;k++)
		for (j=0; j<NumClust; j++) {
			index = j * NumOfVectors/NumClust ; 
			Center[j][k] = (FeatVect+index)->dimen[k];
		}
}



