#include <stdio.h>
#include <stdlib.h>
#include "CA.h"

#include <iostream>
#include <fstream>

int main(int argc, char **argv)
{
	int j,i,k,NumClust, MaximumIt, NumOfVectors, Dim, EXPINC, EXPDEC;
	double Eita_0, *MinPts, *MinPts2, TAU, m_p;
	FILE *fp,*in1;
	char *ca_path;
	char savepath[1024];
	
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
	ca_path      = argv[10];
	
	double **centers;
	double **W;
	int W_i = 0, W_j = 0;

	// Read the input Data
	printf ("Start Reading\n");
	//   	snprintf(savepath,1023,"%s/frames",output_path);
	// in1 = fopen ("/private/tmp/TUF_pjd/super_voxel_model/ca/InData.dat", "rb");
	snprintf( savepath, 1023, "%s/InData.dat", ca_path );
	in1 = fopen( savepath, "rb" );

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
	// in1 = fopen ("/private/tmp/TUF_pjd/super_voxel_model/ca/InCenters.dat", "rb");
	snprintf( savepath, 1023, "%s/InCenters.dat", ca_path );
	in1 = fopen( savepath, "rb" );

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
	// in1 = fopen ("/private/tmp/TUF_pjd/super_voxel_model/ca/MinClSize.dat", "rb");
	snprintf( savepath, 1023, "%s/MinClSize.dat", ca_path );
	printf( "#1:  %s\n", savepath );

	in1 = fopen( savepath, "rb" );
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


	// Read Mahalanobis Distances
	printf ( "Start Reading Feature Co-Variance\n" );
	// in1 = fopen ( "/private/tmp/TUF_pjd/super_voxel_model/ca/MahalDist.dat", "rb" );
	snprintf( savepath, 1023, "%s/MahalDist.dat", ca_path );
	in1 = fopen( savepath, "rb" );
	
	// W is inverse of feature (Dim) co-variance
	W = DMatrix_2D( 0, Dim, 0, Dim );
//	mahals = DMatrix_1D (0, NumOfVectors);	
	double *t3 = DMatrix_1D (0, Dim);
	
	for ( W_i = 0; W_i < Dim; W_i++ ) {
		fread( t3, sizeof(double), Dim, in1 );

		for ( W_j = 0; W_j < Dim; W_j++ ) {
			W[ W_i ][ W_j ] = t3[ W_j ];
		}
	}
	printf("End Reading Feature Co-Variance\n");
	fclose(in1);

	// Call CA algorithm
    NumClust=CA(FeatVect, NumOfVectors, Dim, NumClust, centers, MaximumIt, Eita_0, MinPts, MinPts2, TAU, EXPINC, EXPDEC, m_p, W, ca_path);

	// Write output number of clusters
	// fp=fopen("/private/tmp/TUF_pjd/super_voxel_model/ca/OutClustNum.txt","w");
	snprintf( savepath, 1023, "%s/OutClustNum.txt", ca_path );
	fp = fopen( savepath, "w" );
	fprintf(fp,"%d\n",NumClust);
	fclose(fp);

	// Write output assignment of points
	// fp=fopen("/private/tmp/TUF_pjd/super_voxel_model/ca/OutCluster.txt","w");
	snprintf( savepath, 1023, "%s/OutCluster.txt", ca_path );
	fp = fopen( savepath, "w" );
	for (i=0;i<NumOfVectors;i++)
		fprintf(fp,"%d\n",(FeatVect+i)->cluster);
	fclose(fp);

	// Write output Centers
	// fp=fopen("/private/tmp/TUF_pjd/super_voxel_model/ca/OutCenters.txt","w");
	snprintf( savepath, 1023, "%s/OutCenters.txt", ca_path );
	fp = fopen( savepath, "w" );
	for (j=0;j<NumClust;j++){
		for (i=0;i<Dim;i++){
 			fprintf(fp,"%f ",centers[j][i]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	Free_Feat(FeatVect,NumOfVectors,OrigNumClust,Dim);	
}

