/*
 * Author: Hantao Nie (nht@pku.edu.cn)
 * Date: 2023-08-09 12:01:06
 * LastEditors: Hantao Nie (nht@pku.edu.cn)
 * LastEditTime: 2023-09-11 16:57:23
 * Description: 
 * 
 * Copyright (c) 2023, Hantao Nie, Peking University. 
 */
// compute the differential of a spectral operator
// given a spectral opertor:
// F(Z) = V * diag(phi(d)) * V'  with Z = V * diag(d) * V'
// its differential is computed as 
// DF(Z)[B] = V * (Omega .* (V' * B * V) ) * V'
// where Omega_{ij} = (phi(d_i) - phi(d_j)) / (d_i-d_j)  if d_i ~= d_j
//                    Dphi(d_i)                          if d_i == d_j
// here Dphi is the defferential of phi

// This function is used for constructing the matrix Omega given d, phi(d) and Dphi(d)


#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check number of inputs and outputs
    if (nrhs != 3) {
        mexErrMsgTxt("Three inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgTxt("One output required.");
    }

    // Get the inputs
    double *d = mxGetPr(prhs[0]);
    double *phid = mxGetPr(prhs[1]);
    double *Dphid = mxGetPr(prhs[2]);

    double tol = 1e-16;

    // Get the length of the inputs
    mwSize n = mxGetN(prhs[0]) * mxGetM(prhs[0]);

    // Create the output matrix
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
    double *Omega = mxGetPr(plhs[0]);

    // Compute Omega
    for (mwSize i = 0; i < n; i++) {
        for (mwSize j = 0; j <= i; j++) {
            if (fabs(d[i] - d[j]) > tol) {
                Omega[i + j*n] = (phid[i] - phid[j]) / (d[i] - d[j]);
            } else {
                Omega[i + j*n] = Dphid[i];
            }

            // by symmetry
            if (j != i){
                Omega[j + i*n] = Omega[i + j*n];
            }
        }
    }
}