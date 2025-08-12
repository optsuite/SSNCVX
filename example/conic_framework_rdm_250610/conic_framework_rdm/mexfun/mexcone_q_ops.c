
// % some useful operations for varibles in  Cartesian product of several second order cones
// out = mexcone_q_ops(X, operand, cone_size)
// example:
//   input X with size(X) ==  [14, 1]
//   operand  is chosen from "norm", "det", "sqrtdet", "inv", "Q", "Qinv", "norm_xbar"
//   cone_size = [2 3 4 5]
//   here x is concatenated by four cones
// the size of output is determined by operand
// Based on the given input X and cone_size, the size of the output depends on the chosen operand. Here's a breakdown of the output size for each operand:
// 1. norm: The output will be a column vector with the same number of rows as the number of elements in cone_size. In this case, the size of the output will be [4, 1].
// 2. det: Similar to the "norm" operand, the output will be a column vector with the same number of rows as the number of elements in cone_size. The size of the output will be [4, 1].
// 3. sqrtdet: The output size will be the same as for "det" and "norm" operands, which is [4, 1].
// 4. inv: The output will have the same size as the input X. In this case, the size of the output will be [14, 1].
// 5. Q: This operand is not implemented yet, so the output size cannot be determined.
// 6. Qinv: This operand is not implemented yet, so the output size cannot be determined.
// 7. norm_xbar: The output will be a column vector with the same number of rows as the number of elements in cone_size. In this case, the size of the output will be [4, 1].
// In summary, for the implemented operands, the output size will be either [4, 1] (for "norm", "det", and "sqrtdet") or [14, 1] (for "inv"). 



#include "mex.h"
#include "omp.h"
#include <math.h>
#include "string.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    double * X;
    char * operand;
    double * cone_size_in;
    mwIndex *cone_size;
    mwIndex *cumsum_cone_size, *ind_head;
    mwIndex k; // number of blocks
    mwSize m_in, n_in, n;
    mwIndex i, j;
    double s;


    if (nrhs != 3){
        mexErrMsgTxt("Three inputs needed.");
        return;
    }

    if (nlhs > 1){
        mexErrMsgTxt("Too many output arguments.");
        return;
    }

    // Get X
    m_in = mxGetM(prhs[0]);
    n_in = mxGetN(prhs[0]);
    n = m_in * n_in;
    if ( n_in != 1 )
        mexErrMsgTxt("X must be column vector.");

    X = mxGetDoubles(prhs[0]);
    if (X == NULL || mxIsSparse(prhs[0])){
        mexErrMsgTxt("Argument 1 must be double");
        return;
    }

    // Get operand
    operand = mxArrayToString(prhs[1]);
    if (operand == NULL){
        mexErrMsgTxt("Argument 2 must be string");
        return;
    }

    if (!( strcmp(operand, "norm") == 0|| strcmp(operand, "det") == 0|| strcmp(operand, "sqrtdet") == 0 || strcmp(operand, "inv") == 0|| strcmp(operand, "Q") == 0|| strcmp(operand, "Qinv") == 0 ) || strcmp(operand, "norm_xbar") == 0 )
        mexErrMsgTxt("Argument 2 must be one of 'norm', 'det', 'sqrtdet', 'inv', 'Q', 'Qinv', 'norm_xbar'");

    // Get cone_size
    k = mxGetM(prhs[2]) * mxGetN(prhs[2]);
    cone_size_in = mxGetPr(prhs[2]);

    // transform cone_size to integer array
    cone_size = mxMalloc(k * sizeof(mwIndex));
    #ifndef SERIAL
    #pragma omp parallel for private(j) shared(cone_size, cone_size_in, k)
    #endif
    for (j=0; j<k; j++){
        cone_size[j] = (mwIndex) cone_size_in[j];
    }

    // record the first index of each block
    ind_head = mxMalloc(k * sizeof(mwIndex));
    ind_head[0] = 0;
    for (j=1; j<k; j++){
        ind_head[j] = ind_head[j-1] + cone_size[j-1];
    }
    if ( ind_head[k-1] + cone_size[k-1] != n)
        mexErrMsgTxt("cone_size does not match X");


    if (strcmp(operand, "norm") == 0){
        double s;
        plhs[0] = mxCreateDoubleMatrix(k, 1, mxREAL);
        #ifndef SERIAL
        #pragma omp parallel for private(s, j) shared(cone_size, ind_head, k, X, plhs) 
        #endif
        for (j=0; j<k; j++){
            s = 0;
            for (i=0; i<cone_size[j]; i++){
                s += X[ind_head[j]+i] * X[ind_head[j]+i];
            }
            mxGetPr(plhs[0])[j] = sqrt(s);
        }
    }
    else if (strcmp(operand, "norm_xbar") == 0){
        double s;
        plhs[0] = mxCreateDoubleMatrix(k, 1, mxREAL);
        #ifndef SERIAL
        #pragma omp parallel for private(s, j) shared(cone_size, ind_head, k, X, plhs) 
        #endif
        for (j=0; j<k; j++){
            s = 0;
            for (i=1; i<cone_size[j]; i++){
                s += X[ind_head[j]+i] * X[ind_head[j]+i];
            }
            mxGetPr(plhs[0])[j] = sqrt(s);
        }
    }
    else if (strcmp(operand, "det") == 0) {
        double s;
        plhs[0] = mxCreateDoubleMatrix(k, 1, mxREAL);
        #ifndef SERIAL
        #pragma omp parallel for private(s, j, i) shared(cone_size, ind_head, k, X, plhs) 
        #endif
        for (j=0; j<k; j++){
            s = X[ind_head[j]] * X[ind_head[j]];
            for (i=1; i<cone_size[j]; i++){
                s -= X[ind_head[j]+i] * X[ind_head[j]+i];
            }
            mxGetPr(plhs[0])[j] = s;
        }
    }
    else if (strcmp(operand, "sqrtdet") == 0) {
        double s;
        plhs[0] = mxCreateDoubleMatrix(k, 1, mxREAL);
        #ifndef SERIAL
        #pragma omp parallel for private(s, j, i) shared(cone_size, ind_head, k, X, plhs) 
        #endif
        for (j=0; j<k; j++){
            s = X[ind_head[j]] * X[ind_head[j]];
            for (i=1; i<cone_size[j]; i++){
                s -= X[ind_head[j]+i] * X[ind_head[j]+i];
            }
            mxGetPr(plhs[0])[j] = sqrt(s);
        }
    }
    else if (strcmp(operand, "inv") == 0) {
        double s;
        plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
        #ifndef SERIAL
        #pragma omp parallel for private(s, j, i) shared(cone_size, ind_head, k, X, plhs) 
        #endif
        for (j=0; j<k; j++){
            // compute det
            s = X[ind_head[j]] * X[ind_head[j]];
            for (i=1; i<cone_size[j]; i++){
                s -= X[ind_head[j]+i] * X[ind_head[j]+i];
            }

            if (s <= 0){
                #ifdef __APPLE__
                mexWarnMsgTxt("MEX Warning in computing inv: determinant is not positive.");
                #endif
            }
 
            // compute inverse
            mxGetPr(plhs[0])[ind_head[j]] = X[ind_head[j]] / s;
            for (i=1; i<cone_size[j]; i++){
                mxGetPr(plhs[0])[ind_head[j]+i] = - X[ind_head[j]+i] / s;
            }
        }
    }
    else if (strcmp(operand, "Q") == 0) {
        mexErrMsgTxt("Not implemented yet.");
    }
    else if (strcmp(operand, "Qinv") == 0) {
        mexErrMsgTxt("Not implemented yet.");
    }
    else{
        mexErrMsgTxt("This should not happen.");
    }

    mxFree(cone_size);
    mxFree(ind_head);

}
