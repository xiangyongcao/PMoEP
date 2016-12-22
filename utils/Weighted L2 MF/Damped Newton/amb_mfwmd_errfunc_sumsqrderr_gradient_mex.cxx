/* 
 * amb_mfwmd_errfunc_sumsqrderr_gradient_mex.cxx
 *
 * Calculates the first derivatives of ||W.*(M-A*B')||^2
 *
 * Written by Aeron Buchanan, amb@robots.ox.ac.uk
 */


#include <mex.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

#define for if(0);else for

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
   // function fd = amb_mfwmd_errfunc_sumsqrderr_gradient_mex(A,B,M,W)

   if (nrhs != 4)
      mexErrMsgTxt("amb_mfwmd_errfunc_sumsqrderr_gradient_mex must have 4 arguments");
   if (nlhs != 1)
      mexErrMsgTxt("amb_mfwmd_errfunc_sumsqrderr_gradient_mex must have 1 return value");

   if (!mxIsDouble(prhs[0]))
     mexErrMsgTxt("amb_mfwmd_errfunc_sumsqrderr_gradient_mex takes only double arguments for A");
   if (!mxIsDouble(prhs[1]))
     mexErrMsgTxt("amb_mfwmd_errfunc_sumsqrderr_gradient_mex takes only double arguments for B");
   if (!mxIsDouble(prhs[2]))
     mexErrMsgTxt("amb_mfwmd_errfunc_sumsqrderr_gradient_mex takes only double arguments for M");
   if (!mxIsDouble(prhs[3]))
     mexErrMsgTxt("amb_mfwmd_errfunc_sumsqrderr_gradient_mex takes only double arguments for W");

   mxArray const* A_ptr = prhs[0];
   mxArray const* B_ptr = prhs[1];
   mxArray const* M_ptr = prhs[2];
   mxArray const* W_ptr = prhs[3];

   int m = mxGetM(M_ptr);
   int n = mxGetN(M_ptr);
   int r = mxGetN(A_ptr);

   //printf("amb_wls_gradient_mex: size %d x %d  rank %d\n", m,n,r);
   
   int mr = m*r;
   int nr = n*r;
  
   // make output array
   mxArray* fd_ptr = mxCreateDoubleMatrix(mr+nr, 1, mxREAL);
  
   double const* A_data = mxGetPr(A_ptr);
   double const* B_data = mxGetPr(B_ptr);
   double const* M_data = mxGetPr(M_ptr);
   double const* W_data = mxGetPr(W_ptr);
   
   double* fd_data = mxGetPr(fd_ptr);
   
#define get_from_mx(M,m,n,i,j) ((M)[(j) * (m) + (i)])
#define getA(i,j) get_from_mx(A_data, m, r, i, j)
#define getB(i,j) get_from_mx(B_data, n, r, i, j)
#define getM(i,j) get_from_mx(M_data, m, n, i, j)
#define getW(i,j) get_from_mx(W_data, m, n, i, j)

  
  // compute first half of d
  for(int a = 0; a < m; ++a){
    for(int b = 0; b < r; ++b){
      int index = a*r+b;   
      for(int j = 0; j < n; ++j){
        double AapBjp_sum = 0;
        for(int p=0; p<r; ++p){
          AapBjp_sum += getA(a,p)*getB(j,p);  
        }
        fd_data[index] += getW(a,j)*getW(a,j)*getB(j,b)*(AapBjp_sum - getM(a,j));
      }
      fd_data[index] *= 2;
    }
  }
  // compute second half of d
  for(int c = 0; c < n; ++c){
    for(int d = 0; d < r; ++d){
      int index = mr + c*r+d;
      for(int i = 0; i < m; ++i){
        double AipBcp_sum = 0;
        for(int p=0; p<r; ++p){
          AipBcp_sum += getA(i,p)*getB(c,p);  
        }
        fd_data[index] += getW(i,c)*getW(i,c)*getA(i,d)*(AipBcp_sum - getM(i,c));
      }
      fd_data[index] *= 2;
    }
  }
  
   // assign output array
   plhs[0] = fd_ptr;
}
