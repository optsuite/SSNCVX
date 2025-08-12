/*
 * Author: Hantao Nie (nht@pku.edu.cn)
 * Date: 2023-07-12 15:59:17
 * LastEditors: Hantao Nie (nht@pku.edu.cn)
 * LastEditTime: 2023-07-12 16:36:22
 * Description: 
 * 
 * Copyright (c) 2023, Hantao Nie, Peking University. 
 */

/*
compute the norm of each row and each column of a sparse matrix
[row_norm, col_norm] = row_col_norm(A, rol_norm_opt, col_norm_opt)
Input
    A: sparse matrix of size m x n
    row_norm_opt: 1 for 1-norm, 2 for 2-norm, 3 for inf-norm
    col_norm_opt: 1 for 1-norm, 2 for 2-norm, 3 for inf-norm

Output
    row_norm: norm of each row of A, of size m x 1
    col_norm: norm of each column of A, of size n x 1

*/
#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <omp.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwIndex *Ir, *Jc;
    mwIndex row, col, total=0;
    double *pr, *outMatrixRow, *outMatrixCol, val;
    double normRow, normCol, absval;
    int rowNormOpt, colNormOpt;
    size_t m, n;
    
    /* Check for proper number of arguments. */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt( "MATLAB:mexcpp:nargin",
                "MEX-file requires three input arguments.");
    } else if(nlhs!=2) {
        mexErrMsgIdAndTxt( "MATLAB:mexcpp:nargout",
                "MEX-file requires two output arguments.");
    }
    
    /* Get matrix A */
    pr = mxGetPr(prhs[0]);
    Ir = mxGetIr(prhs[0]);
    Jc = mxGetJc(prhs[0]);
    
    /* Get the size of A */
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    
    /* Get norms options */
    rowNormOpt = (int) mxGetScalar(prhs[1]);
    colNormOpt = (int) mxGetScalar(prhs[2]);

    /* Create a matrix for the return arguments */
    plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
    
    /* Assign pointers to output matrices */
    outMatrixRow = mxGetPr(plhs[0]);
    outMatrixCol = mxGetPr(plhs[1]);

    /* Compute the row norms */
// #pragma omp parallel for private(total, row, val, absval) reduction(+:outMatrixRow[:m], outMatrixCol[:n])
    for (col=0; col<n; col++) {
        for (total=Jc[col]; total<Jc[col+1]; total++) {
            row = Ir[total];
            val = pr[total];
            absval = fabs(val);

            if (rowNormOpt == 1 ){
                outMatrixRow[row] += absval;
            }
            else if (rowNormOpt == 2){
                outMatrixRow[row] +=  val * val;
            } else if (rowNormOpt == 3) {
                outMatrixRow[row] = (outMatrixRow[row] < absval) ? absval : outMatrixRow[row];
            }

            if (colNormOpt == 1 ){
                outMatrixCol[col] += absval;
            }
            else if (colNormOpt == 2){
                outMatrixCol[col] +=  val * val;
            } else if (colNormOpt == 3) {
                outMatrixCol[col] = (outMatrixCol[col] < absval) ? absval : outMatrixCol[col];
            }
        }
    }

    /* If norm is L2, we need to take square root */
    if (rowNormOpt == 2) {
        for (row = 0; row < m; row++) {
            outMatrixRow[row] = sqrt(outMatrixRow[row]);
        }
    }

    if (colNormOpt == 2) {
        for (col = 0; col < n; col++) {
            outMatrixCol[col] = sqrt(outMatrixCol[col]);
        }
    }
}
