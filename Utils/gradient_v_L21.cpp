/*
 *  
 * mex gradient_v_L21.cpp COMPFLAGS="/openmp $COMPFLAGS"
 */
 
#include <limits>
#include <iostream>
#include "mex.h"
#include <vector>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{ 
    if (nrhs != 5)
		mexErrMsgTxt("Usage: [S] = gradient_v_L21(C*A,B,Psi'*S,v,G) where ||CA-Psi'*S*diag(v)*G||_{2,1}.");
    if (nlhs > 1)
		mexErrMsgTxt("Too many output arguments."); 

	// input	
    const double *CA = (double*)mxGetPr(prhs[0]);	
    const double *sums_i = (double*)mxGetPr(prhs[1]);	
    const double *PtS = (double*)mxGetPr(prhs[2]);	
    const double *v = (double*)mxGetPr(prhs[3]);
    const double *G = (double*)mxGetPr(prhs[4]);
      
    //matrix MxN
    const int n = int( mxGetM(prhs[0]));
    const int f = int( mxGetN(prhs[0]));
    const int m = int( mxGetM(prhs[4]));
    
//     mexPrintf("n=%d, f=%d, m=%d\n",n,f,m);
    
    // output
	plhs[0] = mxCreateDoubleMatrix( (mwSize)m, 1, mxREAL);
	double* gv = mxGetPr(plhs[0]);
  
//     std::vector<double> sums_i(f);
    
//     #pragma omp parallel for
//     for(int j=0; j<f; ++j)
//     {
//         double sum_i1 = 0;
//         for(int i=0; i<n; ++i)
//         {
//             sum_i1+=(CA[i + j*n]-B[i + j*n])*(CA[i + j*n]-B[i + j*n]);
//         }
//         sums_i[j] = 1.0/std::sqrt(sum_i1);
//     }
    
    #pragma omp parallel for
    for(int p=0; p<m; ++p)   
    {
//         double th = std::tanh(6*(v[p]-0.5));
//         double der_v = (1-th*th)*3;
            
        double sum_j=0;
        for(int j=0; j<f; ++j)
        {
            double sum_i2 = 0;
            for(int i=0; i<n; ++i)
            {
                double t1 = (CA[i + j*n]/*-B[i + j*n]*/);                
                double t2 = -PtS[i + p*n];
                sum_i2 += t1*t2;
            }
            sum_j += sums_i[j]*sum_i2*G[p + j*m];
        }
        gv[p] = sum_j*v[p];
    }
}
