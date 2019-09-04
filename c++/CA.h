/**************************************************************************************************/
/***  All the variables/prototypes used within the CA algorithm											***/
/**************************************************************************************************/

//#include <Math.h>
//#include <iostream>

//#define MaximumIt 10
//#define EITA_0 5.000000
//#define TAU 10.000000
//#define EXPINC 20.000000
//#define EXPDEC 80.000000

#define MIN_SIMILARITY 0.1

#define INFINITY_1 1.0e35
#define EPS 1.0e-6
// #define MAX_CENT_DIFF 0.01 
#define MAX_CENT_DIFF 0.001

//#define	MinPts 1.0  //For min Fuzzy Card
//#define	MinPts2 2 //For min Crisp Card

#define  MIN(a,b) ((a) < (b)) ? (a) : (b)
#define  MAX(a,b) ((a) > (b)) ? (a) : (b)

static double at, bt, ct;
#define PYTHAG(a, b) ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct = bt/at, at*sqrt (1.0+ct*ct)) : (bt ? (ct = at/bt, bt*sqrt (1.0+ct*ct)):0.0))

#define SIGN(a, b) ((b) >= 0.0 ? fabs (a) : -fabs (a))

struct  Feature_Info {
	double  *memship;
	double  *dist;
	double *dimen;
	unsigned short cluster;
};

int CA (struct Feature_Info *FeatVect, int NumOfVectors, int Dim, int NumClust , double **Center, int MaximumIt, double Eita_0, double *MinPts, double *MinPts2, double TAU, int EXPINC, int EXPDEC, double m_p, double **W, char *ca_path);
void EucDistance (struct Feature_Info *FeatVector, double **Center, int NumClust, int NumOfVectors, int Dim, double m_p);
void AssignPts (struct Feature_Info *FeatVector, int NumClust, int NumOfVectors);
double Get_AggCte (struct Feature_Info *FeatVector, int NumClust, int IterNum, int NumOfVectors, double Eita_0, double TAU, int EXPINC, int EXPDEC);
void CompMem (struct  Feature_Info *FeatVector, int NumClust, double AggCte, int NumOfVectors);
void UpdateNumClusters (struct Feature_Info *FeatVect, int * NumClust, int IterNum, int NumOfVectors, double *MinPts, double *MinPts2);
void FuzzMem (struct Feature_Info *FeatVector, int NumClust, int NumOfVectors);
void FuzzCenters (struct Feature_Info *FeatVect, double ** Center, int NumClust, int NumOfVectors, int Dim);
void InitCenters (double **Center, struct Feature_Info *FeatVect, int NumClust, int NumOfVectors, int Dim);
void VecMatMult ( double *, double **, double *, int, int, int );
double DotProduct (double *, double *, int );

void Error_msg (const char msg[100], int num);
double *DMatrix_1D (int X0, int X1);
void Free_DMatrix_1D (double *Mat, int X0, int X1);
double  **DMatrix_2D (int X0, int X1, int Y0, int Y1);
void Free_DMatrix_2D (double **Mat, int X0, int X1, int Y0, int Y1);
int *IMatrix_1D (int X0, int X1);
void Free_IMatrix_1D (int *Mat, int X0, int X1);
struct Feature_Info *MemFeat (int NumOfVectors, int NC, int Dim);
void Free_Feat (struct Feature_Info  *FeatVect, int NumVect, int NC, int Dim);
