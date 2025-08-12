/*
 * @Author: Hantao Nie (nht@pku.edu.cn)
 * @Date: 2023-05-02 16:43:06
 * LastEditors: Hantao Nie (nht@pku.edu.cn)
 * LastEditTime: 2023-08-02 11:37:56
 * @Description: 
 * 
 * Copyright (c) 2023, Hantao Nie, Peking University. 
 */

/*
 *       Compiler:  mex -O -R2018a spmdiam_impl.c
*/

#include "mex.h"
#include "omp.h"

/**
 * @brief sparse matrix multiply diagonal matrix implementation
 * Y = spmdiag(A, D)
 * A: sparse double
 * D: double vector
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    mwIndex j, i;
    if (nrhs != 2){
        mexErrMsgTxt("Two inputs needed.");
        return;
    }

    if (nlhs > 1){
        mexErrMsgTxt("Too many output arguments.");
        return;
    }

    // parameter checking
    double *A = mxGetDoubles(prhs[0]);
    if (A == NULL || !mxIsSparse(prhs[0])){
        mexErrMsgTxt("Argument 1 must be sparse double.");
        return;
    }

    double *D = mxGetDoubles(prhs[1]);
    if (D == NULL || mxIsSparse(prhs[1])){
        mexErrMsgTxt("Argument 2 must be double");
        return;
    }

    // dimension checking
    mwSize n = mxGetN(prhs[0]);
    mwSize n1 = (mwSize)mxGetNumberOfElements(prhs[1]);

    if (n != n1){
        mexErrMsgTxt("Matrix dimensions must agree.");
        return;
    }

    // perform A * D
    mwIndex *A_ir = mxGetIr(prhs[0]);
    mwIndex *A_jc = mxGetJc(prhs[0]);

    plhs[0] = mxDuplicateArray(prhs[0]);
    double *Y = mxGetDoubles(plhs[0]);

    // loop over columns
    #pragma omp parallel for private(j, i) shared(A_jc, Y, n, D) schedule(dynamic, 1)
    for (j = 0; j < n; ++j){
        for (i = A_jc[j]; i < A_jc[j+1]; ++i)
            Y[i] *= D[j];
    }
}
