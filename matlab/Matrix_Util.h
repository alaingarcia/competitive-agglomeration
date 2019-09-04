/****************************************************************************/
/*   This file contains procedures for manipulating matrices and vectors.   */
/*									    */
/*                                                                          */
/*        		Written By: Frigui, Hichem                          */
/*                            Date: Sept 2nd, 1995                          */
/****************************************************************************/




/*void nrerror( )
 {
   fprintf(stderr,"\nrun-time error...\n");
   fprintf(stderr,"%s\n",error_text);
   fprintf(stderr,"...now exiting to system...\n\n");
   exit(1);
}*/


/****************************************************************************/
/* This procedure performs matrix multiplication: m = a x b                 */
/****************************************************************************/
void MatMult (double **a, double **b, double **m, int dima_1, int dima_2, int dimb_1, int dimb_2 )
{
    int i, j, k;

	for(i=0; i<dima_1; i++)
		for(j=0; j<dimb_2; j++){
		m[i][j] = 0.0;
		for(k=0; k<dima_2; k++)	m[i][j] += a[i][k] * b[k][j];
    }
}

/****************************************************************************/
/* This procedure multiplies a vector by a matrix: m = v x a                */
/****************************************************************************/
void VecMatMult (double *v, double **a, double *m, int dimv, int dima_1, int dima_2)
{
	int i,k;
	for(i=0; i<dima_2; i++){
		m[i] = 0.0;
		for(k=0; k<dimv; k++)	m[i] += v[k] * a[k][i];
	}
}

/****************************************************************************/
/* This procedure multiplies a matrix by a vector: m = a x v                */
/****************************************************************************/
void MatVecMult (double **a, double *v, double *m, int dima_1, int dima_2, int dimv)
{
	int i, k;
	for(i=0; i<dima_1; i++){
		m[i] = 0.0;
		for(k=0; k<dimv; k++)	m[i] += a[i][k] * v[k];
	}
}

/****************************************************************************/
/* This procedure computes the transpose of a matrix.                       */
/****************************************************************************/
void MatTrans (double **a, double **a_t, int dima_1, int dima_2)
{
	int i, j;
	for(i=0; i<dima_2; i++)
		for(j=0; j<dima_1; j++)	a_t[i][j] = a[j][i];
}

/****************************************************************************/
/* This procedure performs matrix subtraction: m = a - b                    */
/****************************************************************************/
void MatSubt(double **a, double **b, double **m, int dima_1, int dima_2, int dimb_1, int dimb_2 )
{
	int i, j;
	for(i=0; i<dima_1; i++)
		for(j=0; j<dimb_2; j++)	m[i][j] = a[i][j] - b[i][j];		
}

/****************************************************************************/
/* This procedure performs matrix addition: m = a + b                       */
/****************************************************************************/
void matrix_add(double **a, double **b, double **m, int dima_1, int dima_2, int dimb_1, int dimb_2 )
{
	int i, j;
	for(i=0; i<dima_1; i++)
		for(j=0; j<dimb_2; j++ )	m[i][j] = a[i][j] + b[i][j];		
}

/****************************************************************************/
/* This function computes the dot product of 2 vectors: v1 and v2.          */
/****************************************************************************/
double DotProduct (double *v1, double *v2, int dim)
{
	int i; double result;
	for (i=0, result=0.0; i<dim; i++) result += v1[i]*v2[i];
	return (result);
}

/****************************************************************************/
/* This procedure computes the outer product of 2 vectors.                  */
/****************************************************************************/
void OuterProd (double *v1, int dim_v1, double *v2, int dim_v2, double **a)
{
	int i, j;
	for (i=0; i<dim_v1; i++)
		for (j=0; j<dim_v2; j++)	a[i][j] = v1[i] * v2[j];
}

/****************************************************************************/
/* This function computes the inverse of a matrix.                          */
/*                                                                          */
/*     INPUT:								    */
/*        a:   A 2-D matrix to be inverted.                                 */
/*        n:   Dimensionality of the matrix a.                              */
/*        l:   Cluster number.						    */
/*        y:   A 3-D matrix that holds the inverse matrices of all clusters.*/
/*                                                                          */
/*     RETURNS:	 							    */
/*        Determinant of the inverted matrix.                               */
/****************************************************************************/
double mat_inv (double **a, double ***y, int n, int l)
{
	double *col, d;
	int i, j, *indx;
	void ludcmp(), lubksb();

	indx=IMatrix_1D(0,N_FEAT);
	col=FMatrix_1D(0,N_FEAT);

	ludcmp(a, n, indx, &d);
	for (j = 0; j < n; j++)
		d *= a[j][j];
	for (j = 0; j < n; j++){
		for (i = 0; i < n; i++)
			col[i] = 0.0;
		col[j] = 1.0;
		lubksb(a, n, indx, col);
		for ( i = 0; i < n; i++)
			y[l][i][j] = col[i];
	}
	return (d);
}
void ludcmp(double **a, int n, int *indx, double *d)
{
	int i, imax, j, k;
	double big, dum, sum, temp, *vv;
 
	vv = FMatrix_1D (0, (n-1));
	*d = 1.0;
	for (i = 0; i < n; i++) {
		big = 0.0;
		for (j = 0; j < n; j++)
			if ((temp = fabs(a[i][j])) > big) big = temp;
		if (fabs(big) <= TINY) nrerror("Singular matrix in routine LUDCMP");
		vv[i] = 1.0/big;
	}
	for (j = 0; j < n; j++){
		for (i = 0; i < j; i++) {
			sum = a[i][j];
			for (k = 0; k < i; k++)
				sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for (i = j; i < n; i++){
			sum = a[i][j];
			for (k = 0; k < j; k++)
				sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
			if ((dum = vv[i]*fabs(sum)) >= big) {
				big = dum;
				imax = i;
 			}
		}
		if (j != imax){
			for (k= 0; k < n; k++){
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			*d = -(*d);
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (fabs(a[j][j]) <= TINY) a[j][j] = TINY;
		if (j != n){
			dum = 1.0 / (a[j][j]);
			for (i = j+1; i < n; i++)
				a[i][j] *= dum;
		}
	}
	Free_DMatrix_1D (vv, 0, n-1);
}

void lubksb(double **a, int n, int *indx, double *b)
{
	int i, ii = 0, ip, j;
	double sum;

	for (i = 0; i < n; i++){
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii)
			for (j = (ii-1); j <= i-1; j++)
				sum -= a[i][j]*b[j];
		else if (sum) ii = (i+1);
			b[i] = sum;
	}
	for (i = (n-1); i >= 0; i--) {
		sum = b[i];
		for (j = i+1; j < n; j++)
			sum -= a[i][j]*b[j];
		b[i] = sum/a[i][i];
	}
}
/****************************************************************************/
/* This procedure find the inverses of all covariance matrices.             */
/*									    */
/*  INPUT:                                                                  */
/*    NumOFClusters: Number of clusters.				    */
/*    Mat_A:       A 3-D array that holds the cov. matrices of all clusters */
/*    det:         A 1-D array that holds the determinant of all cov. mat.  */
/*    Dim:         Dimensionality of the feature space.			    */
/*									    */
/*  OUTPUT:                                                                 */
/*    InvMat_A:    A 3-D array that holds the inverses of all Cov. matrices */
/****************************************************************************/
void InvertMat_DIAG (double ****Mat_A, double ****InvMat_A, double **det, int *DimSets)
{
	int i, j, k, n;
	double tempDet;

	for (n=0; n<NumFeatSets; n++){
		for (i=0; i<NumOfClusters; i++) {
			for (j=0; j<DimSets[n]; j++)
				for (k=0; k<DimSets[n]; k++) InvMat_A[n][i][j][k] = 0.0;

			det[n][i] = 1.0;
			tempDet=1.0;
			for (j=0; j<DimSets[n]; j++){
				tempDet*=(Mat_A[n][i][j][j]+EPS);
//				det[n][i] *= (Mat_A[n][i][j][j]+EPS);
				InvMat_A[n][i][j][j] = 1.0/(Mat_A[n][i][j][j]+EPS);
			}
			det[n][i]=(double)pow((double)tempDet,1.0/(double)DimSets[n]);
		}
	}
}


void InvertMat (double ****Mat_A, double ****InvMat_A, double **det, int *DimSets, int NumClust)
{
	int i,j,k,n,l;
	double **temp_A;  
	double ****buffer;

	buffer=FMatrix_4D(0,NumFeatSets,0,NumClust,0,N_FEAT,0,N_FEAT);
	temp_A=FMatrix_2D(0,N_FEAT,0,N_FEAT);

	for (n=0; n<NumFeatSets; n++){
		for (i=0; i<NumClust; i++) {
			for (j=0; j<DimSets[n]; j++){
				for (k=0; k<DimSets[n]; k++) {
					buffer[n][i][j][k] = 0.0;
				}
			}
		}
	}
	for (n=0; n<NumFeatSets; n++){
		for (i=0; i<NumClust; i++) {
			for (j=0; j<DimSets[n]; j++)
				Mat_A[n][i][j][j]+= EPS;
    
			for (j=0; j<DimSets[n]; j++)
				for (k=0; k<DimSets[n]; k++)	
					temp_A[j][k] = Mat_A[n][i][j][k];
			det[n][i] = mat_inv (temp_A, buffer[n], DimSets[n], i);
			det[n][i]=(double)pow((double)det[n][i],1.0/(double)DimSets[n]);
		}
	}
  
  
	for (n=0; n<NumFeatSets; n++){
		for (i=0; i<NumClust; i++) {
			for (j=0; j<DimSets[n]; j++){
				for (k=0; k<DimSets[n]; k++) {
					InvMat_A[n][i][j][k] = buffer[n][i][j][k];
				}
			}
		}
	}
}		       



/*
void InvertMat2 (Mat_A, InvMat_A, det, NumOfClusters, Dim)
 double ***Mat_A, ***InvMat_A, *det; int NumOfClusters, Dim;
{
   int i, j, k;

  
  for (i=0; i<NumOfClusters; i++) {
     for (j=0; j<Dim; j++)
        for (k=0; k<Dim; k++) InvMat_A[i][j][k] = 0.0;

     det[i] = 1.0;
     for (j=0; j<Dim; j++){
        det[i] *= (Mat_A[i][j][j]+EPS);
        InvMat_A[i][j][j] = 1.0/(Mat_A[i][j][j]+EPS);
     }

   }
} 
*/
/*
void InvertMat1 (Mat_A, InvMat_A, det, NumOfClusters, Dim)
 double ***Mat_A, ***InvMat_A, *det; int NumOfClusters, Dim;
{
  double **temp_A;  int i, j, k;
  temp_A = DMatrix_2D (0, Dim-1, 0, Dim-1);

  for (i=0; i<NumOfClusters; i++) {

    for (j=0; j<Dim; j++)
	if ( fabs(Mat_A[i][j][j]) < EPS){
	    printf ("Mat[%d][%d][%d] = %e\n", i, j, j, Mat_A[i][j][j]);
            Mat_A[i][j][j]=EPS;
	}


     for (j=0; j<Dim; j++)
	for (k=0; k<Dim; k++)	temp_A[j][k] = Mat_A[i][j][k];
     det[i] = mat_inv (temp_A, InvMat_A, Dim, i);

   }

  Free_DMatrix_2D (temp_A, 0, Dim-1, 0, Dim-1);
} 
*/
/****************************************************************************/
/* This procedure find the Eigenvalues and eigenvectors of a symmetric mat. */
/*									    */
/*  INPUT:                                                                  */
/*    a:    symmetrix matrix.				    		    */
/*    n:    Dimensionalty of matrix a					    */
/*									    */
/*  OUTPUT:                                                                 */
/*    v:    Eigenvectors of matrix a					    */
/*    d:    Sorted Eigenvalues of matrix a				    */
/*    nrot: Iteration time.						    */
/****************************************************************************/
#define ROTATE(a, i, j, k, l) g = a[i][j]; h = a[k][l]; a[i][j] = g - s*(h+g*tau); a[k][l] = h + s*(g-h*tau);

void jacobi(double **a, int n, double *d, double **v, int *nrot)
{
	int j, iq, ip, i;
	double **temp_a;
	double tresh, theta, tau, t, sm, s, h, g, c, *b, *z;

	temp_a = FMatrix_2D (0, (n-1), 0, (n-1) );
	b = FMatrix_1D (0, (n-1) );
	z = FMatrix_1D (0, (n-1) );

	for (i=0; i<n; i++)
		for (j=0; j<n; j++)	temp_a[i][j] = a[i][j];
   
	for (ip=0; ip<n; ip++){
		for (iq = 0; iq < n; iq++) v[ip][iq] = 0.0;
		v[ip][ip] = 1.0;
	}
	for (ip=0; ip<n; ip++) {
        b[ip] = d[ip] = a[ip][ip];
        z[ip] = 0.0;
	}
	*nrot = 0;
	for (i = 1; i <= 50; i++) {
		sm = 0.0;
        for (ip = 0; ip < n-1; ip++) {
			for (iq = ip + 1; iq < n; iq++)
				sm += fabs(a[ip][iq]);
        }
        if (sm == 0.0) {
            Free_DMatrix_1D (z, 0, (n - 1) );
            Free_DMatrix_1D (b, 0, (n - 1) );
            eigsrt(d, v, n);
			for ( i = 0; i < n; i++ )
				for ( j = 0; j < n; j++ )
			a[i][j] = temp_a[i][j];
   
			return;
        }
        if (i < 4)
			tresh = 0.2*sm / (n*n);
        else
			tresh = 0.0;
        for (ip = 0; ip < n-1; ip++) {
			for (iq = ip+1; iq < n; iq++) {
				g = 100.0 * fabs(a[ip][iq]);
				if (i > 4 && (double)fabs(d[ip])+g == (double)fabs(d[ip]) && (double)fabs(d[iq])+g == (double)fabs(d[iq]))
					a[ip][iq] = 0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h = d[iq] - d[ip];
					if ((double)fabs(h)+g == (double)fabs(h))
						t = (a[ip][iq]) / h;
					else {
						theta = 0.5*h / (a[ip][iq]);
						t = 1.0 / (fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
                    c = 1.0 / sqrt(1+t*t);
                    s = t * c;
                    tau = s / (1.0+c);
                    h = t * a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq] = 0.0;
                    for (j = 0; j <= ip-1; j++) {
						ROTATE(a, j, ip, j, iq)
                    }
                    for (j = ip+1; j <= iq-1; j++) {
                        ROTATE(a, ip, j, j, iq)
                    }
                    for (j = iq + 1; j < n; j++) {
                        ROTATE(a, ip, j, iq, j)
                    }
                    for (j = 0; j < n; j++) {
                        ROTATE(v, j, ip, j, iq)
                    }
                    ++(*nrot);
				}
			}
        }
        for (ip = 0; ip < n; ip++) {
			b[ip] += z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
        }
	}
	printf("Too many iterations in routine JACOBI");
}

/*******************************
 Eigenvec & Eigenval sorting   *
     input : eig_val & eig_vec *
             dimension         *
    output : eig_val & eig_vec *
*******************************/

int eigsrt(double *d, double **v, int n)
{
	int k, j, i;
	double p;
	for (i = 0; i < n - 1; i++) {
        p = d[k=i];
        for (j = i+1; j < n; j++)
            if (d[j] > p) p = d[k=j];
        if (k != i) {
			d[k] = d[i];
			d[i] = p;
			for (j = 0; j < n; j++) {
				p = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = p;
			}
        }
	}
}
  

