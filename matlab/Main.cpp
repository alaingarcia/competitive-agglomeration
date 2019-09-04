#include <stdio.h>
#include <stdlib.h>
#include "CA.h"

int main(int argc, char **argv)
{
	int j,i,k,NumClust, MaximumIt, NumOfVectors, Dim, EXPINC, EXPDEC;
	double Eita_0, *MinPts, *MinPts2, TAU, m_p;
	FILE *fp,*in1;
	struct Feature_Info *FeatVect;
	NumClust     = atoi(argv[1]);
	NumOfVectors = atoi(argv[2]);
	Dim			 = atoi(argv[3]);
	Eita_0       = atof(argv[4]);
	MaximumIt    = atoi(argv[5]);
	TAU			 = atof(argv[6]);	
	EXPINC		 = atoi(argv[7]);	
	EXPDEC		 = atoi(argv[8]);	
	m_p			 = atoi(argv[9]);	//power for distance
	
	double **centers;

	// Read the input Data
	printf ("Start Reading\n");
	in1 = fopen ("InData.dat", "rb");
	FeatVect=MemFeat(NumOfVectors,NumClust,Dim);
	double *t = DMatrix_1D (0, NumOfVectors);
	for (i = 0; i < Dim; i++) {	
		int nread = fread(t, sizeof(double), NumOfVectors, in1);
		for (j = 0; j < NumOfVectors; j++) {
			(FeatVect+j)->dimen[i]= t[j];
		}
	}
	printf("End Reading\n");
	fclose(in1);

	//fp=fopen("InData.dat","r");
	//FeatVect=MemFeat(NumOfVectors,NumClust,Dim);
	//for (i=0;i<NumOfVectors;i++)
	//	for (k=0;k<Dim;k++)
	//		fscanf(fp,"%lf ",&(FeatVect+i)->dimen[k]);
	//fclose(fp);

	int OrigNumClust = NumClust;

   //InitCenters (centers,FeatVect,NumClust,NumOfVectors,Dim);
	centers = DMatrix_2D(0,NumClust,0,Dim);	
	in1 = fopen ("InCenters.dat", "rb");
	double *t1 = DMatrix_1D (0, NumClust);
	for (k = 0; k < Dim; k++){
		int nread = fread(t1, sizeof(double), NumClust, in1);
		for (j = 0; j < NumClust; j++) {
			centers[j][k] = t1[j];
		}
	}
	printf("End Reading Centers\n");
	fclose(in1);

	MinPts = DMatrix_1D (0, MaximumIt);	
	MinPts2 = DMatrix_1D (0, MaximumIt);
    double *t2 = DMatrix_1D (0, MaximumIt);
	in1 = fopen ("MinClSize.dat", "rb");

	fread(t2, sizeof(double), MaximumIt, in1);

	for (j = 0; j < MaximumIt; j++) {
		MinPts[j] = t2[j];
	}
	fread(t2, sizeof(double), MaximumIt, in1);
	for (j = 0; j < MaximumIt; j++) {
		MinPts2[j] = t2[j];
	}
	printf("End Reading Merging Cond.\n");
	fclose(in1);

	// Call CA algorithm
    NumClust=CA(FeatVect, NumOfVectors, Dim, NumClust, centers, MaximumIt, Eita_0, MinPts, MinPts2, TAU, EXPINC, EXPDEC, m_p);

	// Write output number of clusters
	fp=fopen("OutClustNum.txt","w");
	fprintf(fp,"%d\n",NumClust);
	fclose(fp);

	// Write output assignment of points
	fp=fopen("OutCluster.txt","w");
	for (i=0;i<NumOfVectors;i++)
		fprintf(fp,"%d\n",(FeatVect+i)->cluster);
	fclose(fp);

	// Write output Centers
	fp=fopen("OutCenters.txt","w");
	for (j=0;j<NumClust;j++){
		for (i=0;i<Dim;i++){
 			fprintf(fp,"%f ",centers[j][i]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	Free_Feat(FeatVect,NumOfVectors,OrigNumClust,Dim);	
}

