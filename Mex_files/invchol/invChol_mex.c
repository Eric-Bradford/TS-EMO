/*
 * Description:
 * Computes the inverse of a symmetric matrix using the cholesky
 * decomposition. We exploit the symmetric since P=L*L',
 * inv(P)=inv_L*inv_L', more efficient than inv which uses
 * LU decomposition.
 *
 *
 * Author(s):
 * Eric Blake
 * Black River Systems
 * blake@brsc.com
 * 05/14/2012
 *
 * Usage:
 * B = invChol_mex(A)
 *
 *
 * Notes:
 *
 *
 *
 *
 * See also inv
 */

/* If we are on Linux. */
#if !defined(_WIN32) && !defined(_WIN64)
#define dpotrf dpotrf_
#define spotrf spotrf_
#define dpotri dpotri_
#define spotri spotri_
#endif

#include "mex.h"
#include "blas.h"

/* Function declarations */
void dpotrf( char*, mwSize*, double*, mwSize*, mwSignedIndex* );
void dpotri( char*, mwSize*, double*, mwSize*, mwSignedIndex* );
void spotrf( char*, mwSize*, float*, mwSize*, mwSignedIndex* );
void spotri( char*, mwSize*, float*, mwSize*, mwSignedIndex* );

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    mwIndex i, j; /* indices for populating lower triangle when finished. */
    mwSize n;   /* Size of the matrix. */
    mwSignedIndex info;     /* flag for success in dpotrf, spotrf, dpotri, spotri */
    char *uplo = "U";     /* upper or lower triangle */
    mxClassID type;       /* array type */
    
    /* If we don't pass in the correct number of arguments, throw error. */
    if (nrhs!=1) {
        mexErrMsgIdAndTxt("mexFunction:invChol_mex:numInputs",
                "1 input required: A (square matrix)");
    }
    
    type = mxGetClassID(prhs[0]); /* get the array type */
    
    n = mxGetM(prhs[0]); /* input matrix dimension */
    /* check for symmetric matrix*/
    if (n!=mxGetN(prhs[0])) {
        mexErrMsgIdAndTxt("MATLAB:invChol_mex:matchdims",
                "matrix is not symmetric");
    }
	
	/* check for real matrix */
	if (mxIsComplex(prhs[0])) {
		mexErrMsgIdAndTxt("MATLAB:invChol_mex:iscomplex",
                "matrix must be real");
	}
	
    /* create output matrix (fortran modifies the input so we need to copy the input to output) */
    plhs[0]=mxDuplicateArray(prhs[0]);
    
    /* If we passed in an empty return an empty. */
    if (n==0) {
        return;
    }
    
    /* double precision */
    if(type==mxDOUBLE_CLASS)
    {
        double *B; /* double pointer to input & output matrices*/
        B = mxGetPr(plhs[0]); /* output matrix pointer */
              
        dpotrf( uplo, &n, B, &n, &info ); /* Double Cholesky decomposition */
        
        /* check for success */
        if (info<0) {
            mexErrMsgIdAndTxt("MATLAB:invChol_mex:dpotrf:illegalvalue",
                    "cholesky decomposition failed: illegal value ");
        }
        if (info>0) {
            mexErrMsgIdAndTxt("MATLAB:invChol_mex:dpotrf:notposdef",
                    "cholesky decomposition failed: matrix is not positive definite");
        }
                
        dpotri( uplo, &n, B, &n, &info ); /* Double Inverse using Cholesky decomposition */
        
        /* check for success */
        if (info<0) {
            mexErrMsgIdAndTxt("MATLAB:invChol_mex:dpotri:illegalvalue",
                    "failed to invert: illegal value");
        }
        
        if (info>0) {
            mexErrMsgIdAndTxt("MATLAB:invChol_mex:dpotri:singularmatrix",
                    "failed to invert: a diagonal element was 0");
        }
                
        /* populate the lower triangle */
        for (i=0; i<n; i++) {
            for (j=i+1; j<n; j++) {
                B[n*i+j]=B[j*n+i];
            }
        }
        
    } else if (type==mxSINGLE_CLASS) {
        float *B; /* float pointer to input and output matrices */
        B = (float*) mxGetData(plhs[0]); /* output matrix pointer */
        
        spotrf( uplo, &n, B, &n, &info ); /* Double Cholesky decomposition */
        
        /* check for success */
        if (info<0) {
            mexErrMsgIdAndTxt("MATLAB:invChol_mex:spotrf:illegalvalue",
                    "cholesky decomposition failed: illegal value ");
        }
        if (info>0) {
            mexErrMsgIdAndTxt("MATLAB:invChol_mex:spotrf:notposdef",
                    "cholesky decomposition failed: matrix is not positive definite");
        }
        
        spotri( uplo, &n, B, &n, &info ); /* Double Inverse using Cholesky decomposition */
        
        /* check for success */
        if (info<0) {
            mexErrMsgIdAndTxt("MATLAB:invChol_mex:spotri:illegalvalue",
                    "failed to invert: illegal value");
        }
        
        if (info>0) {
            mexErrMsgIdAndTxt("MATLAB:invChol_mex:spotri:singularmatrix",
                    "failed to invert: a diagonal element was 0");
        }
        
        /* populate the lower triangle*/
        for (i=0; i<n; i++) {
            for (j=i+1; j<n; j++) {
                B[n*i+j]=B[j*n+i];
            }
        }
        
    } else {
        mexErrMsgIdAndTxt("MATLAB:invChol_mex:illegaltype",
                "only single or double matrix inputs are allowed");
    }
    
    
    return;
}