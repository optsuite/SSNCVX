/*
 * @Author: Hantao Nie (nht@pku.edu.cn)
 * @Date: 2023-05-16 16:53:56
 * LastEditors: Hantao Nie (nht@pku.edu.cn)
 * LastEditTime: 2023-11-25 15:06:05
 * @Description: 
 * 
 * Copyright (c) 2023, Hantao Nie, Peking University. 
 */
// the proximal mapping and its derivative wrt z
// usage
//  [X] = mexprox_cone_q(Z, mu, cone_size)
//  [X] = mexprox_cone_q(Z, mu, cone_size, det_reg)
//  [X, DXdiag, DXcoeff, DXrank1] = mexprox_cone_q(Z, mu, cone_size)
//  [X, DXdiag, DXcoeff, DXrank1] = mexprox_cone_q(Z, mu, cone_size, det_reg)

// This function computes the proximal mapping of the negative logarithm of the second-order cone
// and its derivative with respect to z. The input arguments are Z, mu, cone_size, and an optional det_reg.
// The output arguments are X, DXdiag, DXcoeff, and DXrank1.

// output
#include "mex.h"
#include "omp.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    double *Z, *X;
    double *DXcoeff, *DXdiag, *DXrank1;
    double mu;
    double det_reg = 0;
    double * cone_size_in;
    mwIndex *cone_size;
    mwIndex *cumsum_cone_size, *ind_head;
    mwIndex k; // number of blocks
    mwSize m_in, n_in, n;
    mwIndex i, j;

    if (nrhs != 3  && nlhs != 4){
        mexErrMsgTxt("Three or Four inputs needed.");
        return;
    }

    if (nlhs != 1 && nlhs != 4){
        mexErrMsgTxt("Too many output arguments.");
        return;
    }

    // Get Z
    m_in = mxGetM(prhs[0]);
    n_in = mxGetN(prhs[0]);
    n = m_in * n_in;
    if (n_in != 1)
        mexErrMsgTxt("Z must be column vector.");

    Z = mxGetDoubles(prhs[0]);
    if (Z == NULL || mxIsSparse(prhs[0])){
        mexErrMsgTxt("Z must be double");
        return;
    }

    // Get mu
    mu = mxGetScalar(prhs[1]);
    if (mu < 0)
        mexErrMsgTxt("mu must be nonnegative.");

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

    // Get det_reg
    if (nrhs == 4){
        det_reg = mxGetScalar(prhs[3]);
        if (det_reg < 0)
            mexErrMsgTxt("det_reg must be nonnegative.");
    }

    // record the first index of each block
    ind_head = mxMalloc(k * sizeof(mwIndex));
    ind_head[0] = 0;
    for (j=1; j<k; j++){
        ind_head[j] = ind_head[j-1] + cone_size[j-1];
    }
    if ( ind_head[k-1] + cone_size[k-1] != n)
        mexErrMsgTxt("cone_size does not match X");



    // allocate output
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    X = mxGetDoubles(plhs[0]);

    if (nlhs == 4){
        plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
        plhs[2] = mxCreateDoubleMatrix(k, 1, mxREAL);
        plhs[3] = mxCreateDoubleMatrix(n, 1, mxREAL);
        DXdiag = mxGetDoubles(plhs[1]);
        DXcoeff = mxGetDoubles(plhs[2]);
        DXrank1 = mxGetDoubles(plhs[3]);
    } else {
        DXdiag = NULL;
        DXcoeff = NULL;
        DXrank1 = NULL;
    }



    // compute proximal mapping
    double detz, norm2z, norm2zbar, norm2xbar, temp0, Delta, temp1, detx, rho;
    #ifndef SERIAL
    #pragma omp parallel for private(j, i, detz, detx, norm2z, norm2zbar, norm2xbar, temp0, Delta, temp1, rho) shared(cone_size, ind_head, k, Z, X, mu, DXdiag, DXcoeff, DXrank1) 
    #endif
    for (j=0; j<k; j++){
        // // implemention 1
        norm2zbar = 0;
        for (i=1; i<cone_size[j]; i++){
            norm2zbar += Z[ind_head[j]+i] * Z[ind_head[j]+i];
        }
        detz = Z[ind_head[j]] * Z[ind_head[j]] - norm2zbar;
        norm2z = Z[ind_head[j]] * Z[ind_head[j]] + norm2zbar;
        Delta = sqrt( detz * detz + 8 * mu * norm2z + 16 * mu * mu);

        

        if (fabs(Z[ind_head[j]]) < 1e-12){  // Z[ind_head[j]] == 0
            detx = mu;
            rho = 1;
        }
        else{ // Z[ind_head[j]] != 0
            // compute detx and rho = mu / detx
            // temp1 = 0.5 * (detz + Delta);
            if (detz >= 0)   temp1 = 0.5 * (detz + Delta);
            else   temp1 = 4 * mu * (norm2z + 2 * mu) / (Delta - detz); // avoid subtracting two close numbers

            if (fabs(temp1) > 1e-12){ // this implies detx is not too small
                if (Z[ind_head[j]] > 0){
                    detx = 0.5 * (temp1 + sqrt(fmax(temp1 * temp1 - 4 * mu * mu, 0.)));
                    rho = mu / detx;
                }
                else{
                    // detx = 0.5 * (temp1 - sqrt(temp1 * temp1 - 4 * mu * mu));
                    detx = 2 * mu * mu / (temp1 + sqrt(temp1 * temp1 - 4 * mu * mu)); // avoid subtracting two close numbers
                    rho = 0.5 * (temp1 + sqrt(fmax(temp1 * temp1 - 4 * mu * mu, 0.))) / mu;
                }
            }
            else{ // this implies detx is small, i.e., degenerate case
                // in practice we observe that, rho is always near the order 1e-1, no matter how small detx is
                // hence, we compute rho first and then compute detx
                // temp1 = 0.5 * (detz + Delta) / mu;
                if (detz >= 0)   temp1 = 0.5 * (detz + Delta) / mu; 
                else   temp1 = 4  * (norm2z + 2 * mu) / (Delta - detz); // avoid subtracting two close numbers

                if (Z[ind_head[j]] > 0){
                    rho = 2 / (temp1 + sqrt(temp1 * temp1 - 4));
                    detx = mu / rho;
                }
                else{
                    rho = 0.5 * (temp1 + sqrt(temp1 * temp1 - 4));
                    detx = mu / rho;
                }
            }
        }
        mxAssert(detx > 0, "determinant of proximal mapping is not positive");
        if (detx < 1e-16){
            #ifdef __APPLE__
            mexWarnMsgTxt("determinant of proximal mapping is <1e-16");
            #endif
        }
        // impletemention 2: using rho
        for (i=1; i<cone_size[j]; i++){
            X[ind_head[j]+i] = Z[ind_head[j]+i] / ( 1 + rho );
        }

        norm2xbar = 0;
        for (i=1; i<cone_size[j]; i++){
            norm2xbar += X[ind_head[j]+i] * X[ind_head[j]+i];
        }
        X[ind_head[j]] = sqrt(detx + norm2xbar);

        mxAssert(X[ind_head[j]] * X[ind_head[j]] - norm2xbar > 0, "X in not in the interior of the second order cone");


        // modify X[ind_head[j]] to ensure that X is in the interior of the second order cone (not too clse to the margin)
        if (X[ind_head[j]] * X[ind_head[j]] - norm2xbar < 1e-16 && det_reg > 0){
            // mexPrintf("modify Xhead, before: detx = %e, real detx = %e, xhead = %e, norm2xbar = %e\n", detx, X[ind_head[j]] * X[ind_head[j]] - norm2xbar, X[ind_head[j]], norm2xbar);
            X[ind_head[j]] = sqrt(fmax(norm2xbar * (1 + det_reg), norm2xbar + det_reg));
            // mexPrintf("after: real detx = %e, xhead = %e, norm2xbar = %e\n", X[ind_head[j]] * X[ind_head[j]] - norm2xbar, X[ind_head[j]], norm2xbar);
        }

        if (X[ind_head[j]] * X[ind_head[j]] - norm2xbar < 1e-16){
            #ifdef __APPLE__
            mexWarnMsgTxt("computing proximal mapping: real determinant of proximal mapping is <1e-16");
            // mexPrintf("computing proximal mapping: rho = %e, detX = %e, real detX = %e, normX = %e, normZ = %e\n", rho, detx, X[ind_head[j]] * X[ind_head[j]] - norm2xbar, X[ind_head[j]] * X[ind_head[j]] + norm2xbar, norm2z);
            #endif
        }

        if (nlhs == 4){
            if (fabs(rho - 1) < 1e-8){
                if (rho >= 1) rho = 1 + 1e-8;
                else rho = 1 - 1e-8;
            }          

                DXdiag[ind_head[j]] = 1 / (1 - rho);
                DXrank1[ind_head[j]] = X[ind_head[j]] * DXdiag[ind_head[j]];
                for (i=1; i<cone_size[j]; i++){
                    DXdiag[ind_head[j] + i] = 1 / (1 + rho);
                    DXrank1[ind_head[j] + i] = - X[ind_head[j] + i] * DXdiag[ind_head[j] + i];
                }

                DXcoeff[j] = - 2  / ( detx / rho + 2  * ( X[ind_head[j]] * X[ind_head[j]] / (1-rho) + norm2xbar / (1+rho) ) );
            
        }
        
        // // implemention 2 
        // norm2zbar = 0;
        // for (i=1; i<cone_size[j]; i++){
        //     norm2zbar += Z[ind_head[j]+i] * Z[ind_head[j]+i];
        // }
        // detz = Z[ind_head[j]] * Z[ind_head[j]] - norm2zbar;
        // norm2z = Z[ind_head[j]] * Z[ind_head[j]] + norm2zbar;
        // Delta = sqrt( detz * detz + 8 * mu * norm2z + 16 * mu * mu);
        // temp1 = 0.5 * (detz + Delta);
        // if (fabs(Z[ind_head[j]]) < 1e-12){  // Z[ind_head[j]] == 0
        //     detx = mu;
        //     lambda = 1;
        //     X[ind_head[j]] = sqrt(0.25 * norm2zbar + mu);
        //     for (i=1; i<cone_size[j]; i++){
        //         X[ind_head[j]+i] = 0.5 * Z[ind_head[j]+i];
        //     }
        // }
        // else{  
        //     if (Z[ind_head[j]] > 0){
        //         detx = 0.5 * (temp1 + sqrt(temp1 * temp1 - 4 * mu * mu));
        //         lambda = mu / detx;
        //     }
        //     else{
        //         // detx = 0.5 * (temp1 - sqrt(temp1 * temp1 - 4 * mu * mu));
        //         detx = 2 * mu * mu / (temp1 + sqrt(temp1 * temp1 - 4 * mu * mu));
        //         lambda = 0.5 * (temp1 + sqrt(temp1 * temp1 - 4 * mu * mu)) / mu;
        //     }
        //     X[ind_head[j]] = Z[ind_head[j]] / (1 - lambda);
        //     for (i=1; i<cone_size[j]; i++){
        //         X[ind_head[j]+i] = Z[ind_head[j]+i] / (1 + lambda);
        //     }
        // }

        // // implemention 3
        // norm2zbar = 0;
        // for (i=1; i<cone_size[j]; i++){
        //     norm2zbar += (Z[ind_head[j]+i] / Z[ind_head[j]] ) * (Z[ind_head[j]+i] / Z[ind_head[j]] );
        // }
        // detz = 1 - norm2zbar;
        // norm2z = 1 + norm2zbar;
        // Delta = sqrt( detz * detz + 8 * mu / (Z[ind_head[j]] * Z[ind_head[j]]) * norm2z + 16 * ( mu / (Z[ind_head[j]] * Z[ind_head[j]]) ) * ( mu / (Z[ind_head[j]] * Z[ind_head[j]])));
        // temp1 = sqrt( 0.5 * (norm2z + Delta + 4 * mu / (Z[ind_head[j]] * Z[ind_head[j]]) ) );
        // X[ind_head[j]] = 0.5 * Z[ind_head[j]] * ( 1 + temp1);
        // temp0 =  1 + 1 / temp1 ;
        // for (i=1; i<cone_size[j]; i++){
        //     X[ind_head[j]+i] = 0.5 * Z[ind_head[j]+i] * temp0;
        // }
    }


    mxFree(cone_size);
    mxFree(ind_head);

}


