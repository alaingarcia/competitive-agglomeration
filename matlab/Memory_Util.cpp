
/**********************************************************************/
/*  This file contains functions that allocate and free memory.       */
/*                                                                    */
/*         Written by:  Hichem Frigui                                 */
/*               Date:  Sept 3rd, 1995.                               */
/**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "CA.h"
/*-------------------------- Error_msg() ------------------------------*/
/*  This function prints out an error massage and exits to the system. */
/*    Input  :  msg ----- The pointer to the error message.            */
/*              num   ---   Exit number.                               */
/***********************************************************************/
void Error_msg (char msg[100], int num)
{
	printf("%s\n",msg);
	exit (num);
}


/********************************* 1-D MATRIX ********************************/
/*======= INTEGER TYPE ========*/

/*-------------------------- IMatrix_1D() -----------------------*/
/*    This function allocates memory to a 1-D matrix of type     */
/*    Integer,  then return the pointer to the matrix		 */
/*    .		                                     		 */
/*    Input  :  X0 --- The first position of the vector.         */
/*              X1 --- The last position of the vector.          */
/*    Output :  Mat -- The pointer to the vector.               */
/*****************************************************************/
int *IMatrix_1D (int X0, int X1)
{
	int  *Mat;

	Mat = (int *)malloc((unsigned)((X1-X0+1)*(sizeof(int))));;
	if(!Mat) Error_msg ("allocation failure in IMatrix_1D function.", 2);
	return((Mat-=X0));
}

/*----------------------- Free_IMatrix_1D() ---------------------*/
/*    This function frees the memory allocated to 1-D matrix     */
/*    of type Integer.                                           */
/*    Input  :  Mat ----- The pointer to the vector.             */
/*              X0 --- The first position of the vector.         */
/*              X1 --- The last position of the vector.          */
/*****************************************************************/
void Free_IMatrix_1D (int *Mat, int X0, int X1)
{
	free ((char *) (Mat+X0));
}


/*======= CHAR TYPE ========*/

/*------------------------- UCMatrix_1D() -----------------------*/
/*    This function allocates memory to a 1-D matrix of type     */
/*    unsigned char,  then return the pointer to the matrix      */
/*    Input  :  X0 --- The first position of the vector.         */
/*              X1 --- The last position of the vector.          */
/*    Output :  Mat -- The pointer to the vector.                */
/*****************************************************************/
char *CMatrix_1D (int X0, int X1)
{
	char  *Mat;

	Mat = (char *)malloc((unsigned)((X1-X0+1)*(sizeof(char))));
	if(!Mat) Error_msg ("allocation failure in CMatrix_1D function.", 2);
	return((Mat-=X0));
}

/*---------------------- Free_UCMatrix_1D() ---------------------*/
/*    This function frees the memory allocated to 1-D matrix     */
/*    of type unsigned char.                                     */
/*    Input  :  Mat ----- The pointer to the vector.             */
/*              X0 --- The first position of the vector.         */
/*              X1 --- The last position of the vector.          */
/*****************************************************************/
void Free_CMatrix_1D (char *Mat, int X0, int X1)
{
	free ((char *) (Mat+X0));
}

/*======= UNSIGNED CHAR TYPE ========*/

/*------------------------- UCMatrix_1D() -----------------------*/
/*    This function allocates memory to a 1-D matrix of type     */
/*    unsigned char,  then return the pointer to the matrix      */
/*    Input  :  X0 --- The first position of the vector.         */
/*              X1 --- The last position of the vector.          */
/*    Output :  Mat -- The pointer to the vector.                */
/*****************************************************************/
unsigned char *UCMatrix_1D (int X0, int X1)
{
	unsigned char  *Mat;

	Mat = (unsigned char *)malloc((unsigned)((X1-X0+1)*
			  		(sizeof(unsigned char))));;
	if(!Mat) Error_msg ("allocation failure in UCMatrix_1D function.", 2);
	return((Mat-=X0));
}

/*---------------------- Free_UCMatrix_1D() ---------------------*/
/*    This function frees the memory allocated to 1-D matrix     */
/*    of type unsigned char.                                     */
/*    Input  :  Mat ----- The pointer to the vector.             */
/*              X0 --- The first position of the vector.         */
/*              X1 --- The last position of the vector.          */
/*****************************************************************/
void Free_UCMatrix_1D (unsigned char *Mat, int X0, int X1)
{
	free ((char *) (Mat+X0));
}


/*======= DOUBLE TYPE ========*/

/*-------------------------- FMatrix_1D() -----------------------*/
/*    This function allocates memory to a 1-D matrix of type     */
/*    double,  then return the pointer to the matrix		 */
/*    .		                                     		 */
/*    Input  :  X0 --- The first position of the vector.         */
/*              X1 --- The last position of the vector.          */
/*    Output :  Mat -- The pointer to the vector.                */
/*****************************************************************/
double *FMatrix_1D (int X0, int X1)
{
	double  *Mat;

	Mat = (double *)malloc((unsigned)((X1-X0+1)*(sizeof(double))));;
	if(!Mat) Error_msg ("allocation failure in FMatrix_1D function.", 2);
	return((Mat-=X0));
}

/*----------------------- Free_FMatrix_1D() ---------------------*/
/*    This function frees the memory  allocated to 1-D matrix    */
/*    of type double.                                            */
/*    Input  :  Mat ----- The pointer to the vector.             */
/*              X0 --- The first position of the vector.         */
/*              X1 --- The last position of the vector.          */
/*****************************************************************/
void Free_FMatrix_1D (double *Mat, int X0, int X1)
{
	free ((char *) (Mat+X0));
}


/*-------------------------- DMatrix_1D() -----------------------*/
/*    This function allocates memory to a 1-D matrix of type     */
/*    double,  then return the pointer to the matrix		 */
/*    .		                                     		 */
/*    Input  :  X0 --- The first position of the vector.         */
/*              X1 --- The last position of the vector.          */
/*    Output :  Mat -- The pointer to the vector.                */
/*****************************************************************/
double *DMatrix_1D (int X0, int X1)
{
	double  *Mat;

	Mat = (double *)malloc((unsigned)((X1-X0+1)*(sizeof(double))));;
	if(!Mat) Error_msg ("allocation failure in DMatrix_1D function.", 2);
	return((Mat-=X0));
}

/*----------------------- Free_DMatrix_1D() ---------------------*/
/*    This function frees the memory  allocated to 1-D matrix    */
/*    of type double.                                            */
/*    Input  :  Mat ----- The pointer to the vector.             */
/*              X0 --- The first position of the vector.         */
/*              X1 --- The last position of the vector.          */
/*****************************************************************/
void Free_DMatrix_1D (double *Mat, int X0, int X1)
{
	free ((char *) (Mat+X0));
}

/********************************* 2-D MATRIX ********************************/
/*======= INTEGER TYPE ========*/

/*----------------------- IMatrix_2D() -------------------------*/
/*   This function allocates memory to a 2-D matrix of type     */
/*   Integer, then returns the pointer to the matrix.           */
/*    Input  :  X0 --- 	The first row of the matrix.            */
/*              X1 --- 	The last row of the matrix.             */
/*              Y0 --- 	The first column of the matrix.         */
/*              Y1 --- 	The last column of the matrix.          */
/*    Output :  Mat --- The pointer to the matrix.              */
/****************************************************************/
int  **IMatrix_2D (int X0, int X1, int Y0, int Y1)
{
	int x, **Mat, *IMtarix_1D();

	Mat = (int **)malloc((unsigned)((X1-X0+1)*(sizeof(int *))));
	if(!Mat) Error_msg ("allocation failure in IMatrix_2D function.",2);
	Mat -= X0;
	for (x=X0; x<=X1; x++) Mat[x] = IMatrix_1D (Y0,Y1);
	return (Mat);
}

/*--------------------- Free_IMatrix_2D() -----------------------*/
/*    This function frees the memory allocated to a 2-D matrix   */
/*    of type integer.                                           */
/*    Input  :  Mat ----- The pointer to the matrix.             */
/*              X0 --- The first row of the matrix.              */
/*              X1 --- The last row of the matrix.               */
/*              Y0 --- The first column of the matrix.           */
/*              Y1 --- The last column of the matrix.            */
/*****************************************************************/
void Free_IMatrix_2D (int **Mat, int X0, int X1, int Y0, int Y1)
{
	int x;

	for (x=X1;x>=X0;x--) free((char *) (Mat[x]+Y0));
	free ((char *)(Mat+X0));
}


/*======= CHAR TYPE ========*/

/*----------------------- UCMatrix_2D() ------------------------*/
/*   This function allocates memory to a 2-D matrix of type     */
/*   unsigned char, then returns the pointer to the matrix.     */
/*    Input  :  X0 --- 	The first row of the matrix.            */
/*              X1 --- 	The last row of the matrix.             */
/*              Y0 --- 	The first column of the matrix.         */
/*              Y1 --- 	The last column of the matrix.          */
/*    Output :  Mat --- The pointer to the matrix.              */
/****************************************************************/
char  **CMatrix_2D (int X0, int X1, int Y0, int Y1)
{
	char  **Mat, *CMtarix_1D(int,int);
	int x;

	Mat = (char **)malloc((unsigned)((X1-X0+1)*
									(sizeof( char *))));
	if(!Mat) Error_msg ("allocation failure in CMatrix_2D function.",2);
	Mat -= X0;
	for (x=X0; x<=X1; x++) Mat[x] = CMatrix_1D (Y0,Y1);
	return (Mat);
}

/*--------------------- Free_UCMatrix_2D() ----------------------*/
/*    This function frees the memory allocated to a 2-D matrix   */
/*    of type unsigned char.                                     */
/*    Input  :  Mat ----- The pointer to the matrix.             */
/*              X0 --- The first row of the matrix.              */
/*              X1 --- The last row of the matrix.               */
/*              Y0 --- The first column of the matrix.           */
/*              Y1 --- The last column of the matrix.            */
/*****************************************************************/
void Free_CMatrix_2D (char **Mat, int X0, int X1, int Y0, int Y1)
{
	int x;

	for (x=X1;x>=X0;x--) free((char *) (Mat[x]+Y0));
	free ((char *)(Mat+X0));
}


/*======= UNSIGNED CHAR TYPE ========*/

/*----------------------- UCMatrix_2D() ------------------------*/
/*   This function allocates memory to a 2-D matrix of type     */
/*   unsigned char, then returns the pointer to the matrix.     */
/*    Input  :  X0 --- 	The first row of the matrix.            */
/*              X1 --- 	The last row of the matrix.             */
/*              Y0 --- 	The first column of the matrix.         */
/*              Y1 --- 	The last column of the matrix.          */
/*    Output :  Mat --- The pointer to the matrix.              */
/****************************************************************/
unsigned char  **UCMatrix_2D (int X0, int X1, int Y0, int Y1)
{
	unsigned char  **Mat, *UCMtarix_1D(int,int);
	int x;

	Mat = (unsigned char **)malloc((unsigned)((X1-X0+1)*
									(sizeof(unsigned char *))));
	if(!Mat) Error_msg ("allocation failure in UCMatrix_2D function.",2);
	Mat -= X0;
	for (x=X0; x<=X1; x++) Mat[x] = UCMatrix_1D (Y0,Y1);
	return (Mat);
}

/*--------------------- Free_UCMatrix_2D() ----------------------*/
/*    This function frees the memory allocated to a 2-D matrix   */
/*    of type unsigned char.                                     */
/*    Input  :  Mat ----- The pointer to the matrix.             */
/*              X0 --- The first row of the matrix.              */
/*              X1 --- The last row of the matrix.               */
/*              Y0 --- The first column of the matrix.           */
/*              Y1 --- The last column of the matrix.            */
/*****************************************************************/
void Free_UCMatrix_2D (unsigned char **Mat, int X0, int X1, int Y0, int Y1)
{
	int x;

	for (x=X1;x>=X0;x--) free((char *) (Mat[x]+Y0));
	free ((char *)(Mat+X0));
}


/*======= double TYPE ========*/

/*----------------------- FMatrix_2D() -------------------------*/
/*   This function allocates memory to a 2-D matrix of type     */
/*   double, then returns the pointer to the matrix.            */
/*    Input  :  X0 --- 	The first row of the matrix.            */
/*              X1 --- 	The last row of the matrix.             */
/*              Y0 --- 	The first column of the matrix.         */
/*              Y1 --- 	The last column of the matrix.          */
/*    Output :  Mat --- The pointer to the matrix.              */
/****************************************************************/
double  **FMatrix_2D (int X0, int X1, int Y0, int Y1)
{
	double **Mat, *FMtarix_1D(int,int);
	int x;

	Mat = (double **)malloc((unsigned)((X1-X0+1)*(sizeof(double *))));
	if(!Mat) Error_msg ("allocation failure in FMatrix_2D function.",2);
	Mat -= X0;
	for (x=X0; x<=X1; x++) Mat[x] = FMatrix_1D (Y0,Y1);
	return (Mat);
}

/*--------------------- Free_FMatrix_2D() -----------------------*/
/*    This function frees the memory allocated to a 2-D matrix   */
/*    of type double.                                            */
/*    Input  :  Mat ----- The pointer to the matrix.             */
/*              X0 --- The first row of the matrix.              */
/*              X1 --- The last row of the matrix.               */
/*              Y0 --- The first column of the matrix.           */
/*              Y1 --- The last column of the matrix.            */
/*****************************************************************/
void Free_FMatrix_2D (double **Mat, int X0, int X1, int Y0, int Y1)
{
	int x;

	for (x=X1;x>=X0;x--) free((char *) (Mat[x]+Y0));
	free ((char *)(Mat+X0));
}

/*======= DOUBLE TYPE ========*/

/*----------------------- DMatrix_2D() -------------------------*/
/*   This function allocates memory to a 2-D matrix of type     */
/*   double, then returns the pointer to the matrix.            */
/*    Input  :  X0 --- 	The first row of the matrix.            */
/*              X1 --- 	The last row of the matrix.             */
/*              Y0 --- 	The first column of the matrix.         */
/*              Y1 --- 	The last column of the matrix.          */
/*    Output :  Mat --- The pointer to the matrix.              */
/****************************************************************/
double  **DMatrix_2D (int X0, int X1, int Y0, int Y1)
{
	double **Mat, *DMtarix_1D(int,int);
	int x;

	Mat = (double **)malloc((unsigned)((X1-X0+1)*(sizeof(double *))));
	if(!Mat) Error_msg ("allocation failure in DMatrix_2D function.",2);
	Mat -= X0;
	for (x=X0; x<=X1; x++) Mat[x] = DMatrix_1D (Y0,Y1);
	return (Mat);
}

/*--------------------- Free_DMatrix_2D() -----------------------*/
/*    This function frees the memory allocated to a 2-D matrix   */
/*    of type double.                                            */
/*    Input  :  Mat ----- The pointer to the matrix.             */
/*              X0 --- The first row of the matrix.              */
/*              X1 --- The last row of the matrix.               */
/*              Y0 --- The first column of the matrix.           */
/*              Y1 --- The last column of the matrix.            */
/*****************************************************************/
void Free_DMatrix_2D (double **Mat, int X0, int X1, int Y0, int Y1)
{
	int x;

	for (x=X1;x>=X0;x--) free((char *) (Mat[x]+Y0));
	free ((char *)(Mat+X0));
}

/********************************* 3-D MATRIX ********************************/
int  ***IMatrix_3D (int X0, int X1, int Y0, int Y1, int Z0, int Z1)
{
	int ***Mat, **IMtarix_2D(int,int,int,int);
	int x;

	Mat = (int ***)malloc((unsigned)((X1-X0+1)*(sizeof(int **))));
	if(!Mat) Error_msg ("allocation failure in IMatrix_3D function.",2);
	Mat -= X0;
	for (x=X0; x<=X1; x++) Mat[x] = IMatrix_2D (Y0, Y1, Z0, Z1);
	return (Mat);
}

/*--------------------- Free_DMatrix_3D() -----------------------*/
/*    This function frees the memory allocated to a 3-D matrix   */
/*    of type double.                                            */
/*    Input  :  Mat ----- The pointer to the matrix.             */
/*              X0 --- The first row of the matrix.              */
/*              X1 --- The last row of the matrix.               */
/*              Y0 --- The first column of the matrix.           */
/*              Y1 --- The last column of the matrix.            */
/*              Z0 --- The first number of the third dimension   */
/*              Z1 --- The last number of the first dimension    */
/*****************************************************************/
void Free_IMatrix_3D (int ***Mat, int X0, int X1, int Y0, int Y1, int Z0, int Z1)
{
	int x1, x2;

	for (x1=X1;x1>=X0;x1--){
		for(x2=Y1;x2>=Y0;x2--)   
		free ((char *) (Mat[x1][x2]+Z0));
		free ((char *) (Mat[x1]+Y0));
	}
	free ((char *)(Mat+X0));
}


/*======= double TYPE ========*/

/*----------------------- FMatrix_3D() -------------------------*/
/*   This function allocates memory to a 3-D matrix of type     */
/*   double, then returns the pointer to the matrix.            */
/*    Input  :  X0 --- 	The first row of the matrix.            */
/*              X1 --- 	The last row of the matrix.             */
/*              Y0 --- 	The first column of the matrix.         */
/*              Y1 --- 	The last column of the matrix.          */
/*              Z0 ---  The first number of the third dimension */
/*              Z1 ---  The last number of the first dimension  */
/*    Output :  Mat --- The pointer to the matrix.              */
/****************************************************************/
double  ***FMatrix_3D (int X0, int X1, int Y0, int Y1, int Z0, int Z1)
{
	double ***Mat, **DMtarix_2D(int,int,int,int);
	int x;

	Mat = (double ***)malloc((unsigned)((X1-X0+1)*(sizeof(double **))));
	if(!Mat) Error_msg ("allocation failure in FMatrix_3D function.",2);
	Mat -= X0;
	for (x=X0; x<=X1; x++) Mat[x] = FMatrix_2D (Y0, Y1, Z0, Z1);
	return (Mat);
}

/*--------------------- Free_FMatrix_3D() -----------------------*/
/*    This function frees the memory allocated to a 3-D matrix   */
/*    of type double.                                            */
/*    Input  :  Mat ----- The pointer to the matrix.             */
/*              X0 --- The first row of the matrix.              */
/*              X1 --- The last row of the matrix.               */
/*              Y0 --- The first column of the matrix.           */
/*              Y1 --- The last column of the matrix.            */
/*              Z0 --- The first number of the third dimension   */
/*              Z1 --- The last number of the first dimension    */
/*****************************************************************/
void Free_FMatrix_3D (double ***Mat, int X0, int X1, int Y0, int Y1, int Z0, int Z1)
{
	int x1, x2;

	for (x1=X1;x1>=X0;x1--){
		for(x2=Y1;x2>=Y0;x2--)   
		free ((char *) (Mat[x1][x2]+Z0));
		free ((char *) (Mat[x1]+Y0));
	}
	free ((char *)(Mat+X0));
}

/*======= DOUBLE TYPE ========*/

/*----------------------- DMatrix_3D() -------------------------*/
/*   This function allocates memory to a 3-D matrix of type     */
/*   double, then returns the pointer to the matrix.            */
/*    Input  :  X0 --- 	The first row of the matrix.            */
/*              X1 --- 	The last row of the matrix.             */
/*              Y0 --- 	The first column of the matrix.         */
/*              Y1 --- 	The last column of the matrix.          */
/*              Z0 ---  The first number of the third dimension */
/*              Z1 ---  The last number of the first dimension  */
/*    Output :  Mat --- The pointer to the matrix.              */
/****************************************************************/
double  ***DMatrix_3D (int X0, int X1, int Y0, int Y1, int Z0, int Z1)
{
	double ***Mat, **DMatrix_2D(int,int,int,int);
	int x;

	Mat = (double ***)malloc((unsigned)((X1-X0+1)*(sizeof(double **))));
	if(!Mat) Error_msg ("allocation failure in DMatrix_3D function.",2);
	Mat -= X0;
	for (x=X0; x<=X1; x++) Mat[x] = DMatrix_2D (Y0, Y1, Z0, Z1);
	return (Mat);
}

/*--------------------- Free_DMatrix_3D() -----------------------*/
/*    This function frees the memory allocated to a 3-D matrix   */
/*    of type double.                                            */
/*    Input  :  Mat ----- The pointer to the matrix.             */
/*              X0 --- The first row of the matrix.              */
/*              X1 --- The last row of the matrix.               */
/*              Y0 --- The first column of the matrix.           */
/*              Y1 --- The last column of the matrix.            */
/*              Z0 --- The first number of the third dimension   */
/*              Z1 --- The last number of the first dimension    */
/*****************************************************************/
void Free_DMatrix_3D (double ***Mat, int X0, int X1, int Y0, int Y1, int Z0, int Z1)
{
	int x1, x2;

	for (x1=X1;x1>=X0;x1--){
		for(x2=Y1;x2>=Y0;x2--)   
		free ((char *) (Mat[x1][x2]+Z0));
		free ((char *) (Mat[x1]+Y0));
	}
	free ((char *)(Mat+X0));
}


/****************************** FEATURE VECTORS ******************************/

/*--------------------------- MemFeat() -------------------------*/
/*    This function allocates memory to the feature vector. The  */
/*    structure of the feature vector is defined in "cluster.h"  */
/*    Input   : NumOfVectors -- Number of feature vectors.       */
/*              NC           -- Number of clusters.		 */
/*    Output  : FeatVect ----   Pointer to the feature vector.   */
/*****************************************************************/

struct Feature_Info *MemFeat (int NumOfVectors, int NC, int Dim)
{
	struct Feature_Info  *Features;
	int i;
  
	Features = (struct Feature_Info *)malloc( (NumOfVectors) * 
					    sizeof(struct Feature_Info ) );
	if(!Features) Error_msg ("allocation failure for feature vectors.",2);
   
    for (i=0; i<NumOfVectors; i++){
		Features[i].memship =  DMatrix_1D (0, NC);
		Features[i].dist    =  DMatrix_1D (0, NC);
		Features[i].dimen   =  DMatrix_1D(0, Dim);
    }

	return (Features);
} 

/*--------------------------- Free_Feat() -----------------------*/
/*    This function frees the memory allocated to the feature    */
/*    vector.                                                    */
/*    Input  :  FeatVect ---- The pointer to the feature vector. */
/*              NC         -- Number of clusters.		 */
/*****************************************************************/

void Free_Feat (struct Feature_Info  *FeatVect, int NumVect, int NC, int Dim)
{
	int i;
	 
	for (i=0; i<NumVect; i++){
		Free_DMatrix_1D (FeatVect[i].dimen,0,Dim);
		Free_DMatrix_1D (FeatVect[i].memship, 0, NC);
		Free_DMatrix_1D (FeatVect[i].dist, 0, NC);
	}  
	free (  (FeatVect));
} 





unsigned char ***UCMatrix_3D(int x0, int x1, int y0, int y1, int z0, int z1)
{
	unsigned char ***Mat, **UCMatrix_2D(int,int,int,int);
	int x;

	Mat=(unsigned char ***)malloc((unsigned)((x1-x0+1)*(sizeof(unsigned char **))));
	if(!Mat) Error_msg ("Allocation failure in UCMatrix_3D function.",2);
	Mat-=x0;
	for (x=x0;x<=x1;x++) Mat[x]=UCMatrix_2D(y0,y1,z0,z1);
	return(Mat);
}




/*----------------------- DMatrix_4D() -------------------------*/
/*   This function allocates memory to a 4-D matrix of type     */
/*   double, then returns the pointer to the matrix.            */
/*    Input  :  X0 --- 	The first row of the matrix.            */
/*              X1 --- 	The last row of the matrix.             */
/*              Y0 --- 	The first column of the matrix.         */
/*              Y1 --- 	The last column of the matrix.          */
/*              Z0 ---  The first number of the third dimension */
/*              Z1 ---  The last number of the first dimension  */
/*              R0 ---  The first number of the fourth dimension*/
/*              R1 ---  The last number in the fourth dimension */
/*    Output :  Mat --- The pointer to the matrix.              */
/****************************************************************/
double  ****DMatrix_4D (int X0, int X1, int Y0, int Y1, int Z0, int Z1, int R0, int R1)
{
	double ****Mat, ***DMatrix_3D(int,int,int,int,int,int);
	int x;

	Mat = (double ****)malloc((unsigned)((X1-X0+1)*(sizeof(double ***))));
	if(!Mat) Error_msg ("allocation failure in DMatrix_4D function.",2);
	Mat -= X0;
	for (x=X0; x<=X1; x++) Mat[x] = DMatrix_3D (Y0, Y1, Z0, Z1,R0 ,R1 );
	return (Mat);
}

/*--------------------- Free_DMatrix_4D() -----------------------*/
/*    This function frees the memory allocated to a 3-D matrix   */
/*    of type double.                                            */
/*    Input  :  Mat ----- The pointer to the matrix.             */
/*              X0 --- The first row of the matrix.              */
/*              X1 --- The last row of the matrix.               */
/*              Y0 --- The first column of the matrix.           */
/*              Y1 --- The last column of the matrix.            */
/*              Z0 --- The first number of the third dimension   */
/*              Z1 --- The last number of the first dimension    */
/*****************************************************************/
void Free_DMatrix_4D (double ****Mat, int X0, int X1, int Y0, int Y1, int Z0, int Z1, int R0, int R1)
{
	int x1, x2, x3;
	  
	for (x1=X1;x1>=X0;x1--){
		for(x2=Y1;x2>=Y0;x2--){
		for (x3=Z1;x3>=Z0;x3--)
		free ((char *) (Mat[x1][x2][x3]+R0));
		free ((char *) (Mat[x1][x2]+Z0));
		}
		free ((char *) (Mat[x1]+Y0));
	}
	free ((char *)(Mat+X0));
}





/*----------------------- FMatrix_4D() -------------------------*/
/*   This function allocates memory to a 4-D matrix of type     */
/*   double , then returns the pointer to the matrix.            */
/*    Input  :  X0 --- 	The first row of the matrix.            */
/*              X1 --- 	The last row of the matrix.             */
/*              Y0 --- 	The first column of the matrix.         */
/*              Y1 --- 	The last column of the matrix.          */
/*              Z0 ---  The first number of the third dimension */
/*              Z1 ---  The last number of the first dimension  */
/*              R0 ---  The first number of the fourth dimension*/
/*              R1 ---  The last number in the fourth dimension */
/*    Output :  Mat --- The pointer to the matrix.              */
/****************************************************************/
double  ****FMatrix_4D (int X0, int X1, int Y0, int Y1, int Z0, int Z1, int R0, int R1)
{
	double ****Mat, ***FMatrix_3D(int,int,int,int,int,int);
	int x;

	Mat = (double ****)malloc((unsigned)((X1-X0+1)*(sizeof(double ***))));
	if(!Mat) Error_msg ("allocation failure in FMatrix_4D function.",2);
	Mat -= X0;
	for (x=X0; x<=X1; x++) Mat[x] = FMatrix_3D (Y0, Y1, Z0, Z1,R0 ,R1 );
	return (Mat);
}

/*--------------------- Free_FMatrix_4D() -----------------------*/
/*    This function frees the memory allocated to a 3-D matrix   */
/*    of type double.                                            */
/*    Input  :  Mat ----- The pointer to the matrix.             */
/*              X0 --- The first row of the matrix.              */
/*              X1 --- The last row of the matrix.               */
/*              Y0 --- The first column of the matrix.           */
/*              Y1 --- The last column of the matrix.            */
/*              Z0 --- The first number of the third dimension   */
/*              Z1 --- The last number of the first dimension    */
/*****************************************************************/
void Free_FMatrix_4D (double ****Mat, int X0, int X1, int Y0, int Y1, int Z0, int Z1, int R0, int R1)
{
	int x1, x2, x3;
	  
	for (x1=X1;x1>=X0;x1--){
		for(x2=Y1;x2>=Y0;x2--){
		for (x3=Z1;x3>=Z0;x3--)
		free ((char *) (Mat[x1][x2][x3]+R0));
		free ((char *) (Mat[x1][x2]+Z0));
		}
		free ((char *) (Mat[x1]+Y0));
	}
	free ((char *)(Mat+X0));
}
