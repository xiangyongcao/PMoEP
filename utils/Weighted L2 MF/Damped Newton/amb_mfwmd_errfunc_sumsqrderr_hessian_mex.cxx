/* 
 * amb_mfwmd_errfunc_sumsqrderr_hessian_mex.cxx
 *
 * Calculates the second derivatives of ||W.*(M-A*B')||^2
 *
 * Written by Aeron Buchanan, amb@robots.ox.ac.uk
 */

#include <mex.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

#define for if(0);else for

// #define amb_debug
#ifndef amb_debug
   #define get_from_mx(M,m,n,i,j) ((M)[(j) * (m) + (i)])
#else
// this will require slightly more thought as 'get_from_mx' is used as an l-value
inline double get_from_mx(double const * M, int m, int n, int i, int j)
{
   if (i < 0 || j < 0 || i >= m || j >= n) {
     printf("error: attempt to access element (%d,%d) of %dx%d matrix\n", i,j,m,n);
     mexErrMsgTxt("bounds error");
   }
   return ((M)[(j) * (m) + (i)]);
}
#endif

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
   if (nrhs != 4)
      mexErrMsgTxt("amb_mfwmd_errfunc_sumsqrderr_hessian_mex must have 4 arguments");
   if (nlhs != 3)
      mexErrMsgTxt("amb_mfwmd_errfunc_sumsqrderr_hessian_mex must have 3 return values");

   if (!mxIsDouble(prhs[0]))
     mexErrMsgTxt("amb_mfwmd_errfunc_sumsqrderr_hessian_mex takes only double arguments for A");
   if (!mxIsDouble(prhs[1]))
     mexErrMsgTxt("amb_mfwmd_errfunc_sumsqrderr_hessian_mex takes only double arguments for B");
   if (!mxIsDouble(prhs[2]))
     mexErrMsgTxt("amb_mfwmd_errfunc_sumsqrderr_hessian_mex takes only double arguments for M");
   if (!mxIsDouble(prhs[3]))
     mexErrMsgTxt("amb_mfwmd_errfunc_sumsqrderr_hessian_mex takes only double arguments for W");

   mxArray const* A_ptr = prhs[0];
   mxArray const* B_ptr = prhs[1];
   mxArray const* M_ptr = prhs[2];
   mxArray const* W_ptr = prhs[3];

   int m = mxGetM(M_ptr);
   int n = mxGetN(M_ptr);
   int r = mxGetN(A_ptr);

   int mr = m*r;
   int nr = n*r;
  
   // make output array
   mxArray* HA_ptr = mxCreateDoubleMatrix(mr, mr, mxREAL);
   mxArray* HB_ptr = mxCreateDoubleMatrix(mr, nr, mxREAL);
   mxArray* HC_ptr = mxCreateDoubleMatrix(nr, r, mxREAL);
  
   double const* A_data = mxGetPr(A_ptr);
   double const* B_data = mxGetPr(B_ptr);
   double const* M_data = mxGetPr(M_ptr);
   double const* W_data = mxGetPr(W_ptr);
   
   double* HA_data = mxGetPr(HA_ptr);
   double* HB_data = mxGetPr(HB_ptr);
   double* HC_data = mxGetPr(HC_ptr);
   
#define getA(i,j) get_from_mx(A_data, m, r, i, j)
#define getB(i,j) get_from_mx(B_data, n, r, i, j)
#define getM(i,j) get_from_mx(M_data, m, n, i, j)
#define getW(i,j) get_from_mx(W_data, m, n, i, j)
#define getHA(i,j) get_from_mx(HA_data, mr, mr, i, j)
#define getHB(i,j) get_from_mx(HB_data, mr, nr, i, j)
#define getHC(i,j) get_from_mx(HC_data, nr, r, i, j)

   /// compute upper-left block of H: HA
   for(int a = 0; a < m; ++a){
     for(int b = 0; b < r; ++b){
       int second_index = a*r+b;
       for(int f = 0; f < r; ++f){
         int first_index = a*r+f;
         for(int j = 0; j < n; ++j){ 
           getHA(first_index,second_index) += getW(a,j)*getW(a,j)*getB(j,b)*getB(j,f);
         }
         getHA(first_index,second_index) *= 2;
       }
     }
   }
   
   // compute upper-right block of H: HB
   for(int c = 0; c < n; ++c){
     for(int d = 0; d < r; ++d){
       //printf("at (%d,%d)\n", c,d);
       int second_index = c*r+d;
       for(int e = 0; e < m; ++e){
         for(int f = 0; f < r; ++f){
           int first_index = e*r+f;
           double AepBcp_sum = 0;
           if (d==f) {
             for(int p = 0; p < r; ++p){
               AepBcp_sum += getA(e,p)*getB(c,p);  
             }
             AepBcp_sum = AepBcp_sum - getM(e,c);
           }
           getHB(first_index,second_index) = 2*getW(e,c)*getW(e,c)*(getA(e,d)*getB(c,f) + AepBcp_sum); 
         }
       }
     }
   }
   
   // printf("done hb\n");
   // compute lower-right block of H: HC
   for(int c = 0; c < n; ++c){
     for(int d = 0; d < r; ++d){
       int second_index = d; // slide diagonal blocks of HC to first r columns
       for(int h = 0; h < r; ++h){
         int first_index = c*r+h;
         for(int i = 0; i < m; ++i){
           getHC(first_index,second_index) += getW(i,c)*getW(i,c)*getA(i,d)*getA(i,h);
         }
         getHC(first_index,second_index) *= 2;
       }
     }
   }
   
   // assign output array
   plhs[0] = HA_ptr;
   plhs[1] = HB_ptr;
   plhs[2] = HC_ptr;
}
