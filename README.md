README
================================================================================================
Version 1.0, 29-Dec-2015

This package contains the MATLAB implementation of "Low-Rank Matrix Factorization Under General Mixture Noise Distributions".

The code has been tested with MATLAB 2014b on a PC with 64-bit windows 7.

================================================================================================

Use of this code is free for research purposes only.

================================================================================================

Reference:

Xiangyong Cao, Yang Chen, Qian Zhao, Deyu Meng, Yao Wang, Dong Wang and Zongben Xu,
Low-Rank Matrix Factorization Under General Mixture Noise Distributions, 
15th International Conference on Computer Vision (ICCV), Chile, Dec. 2015 (Oral)

================================================================================================

Installation:

1. Unpack the contents of the compressed file to a new directory.

2. Run the Demos

================================================================================================

Demos:

Demo_EP.m          % EP 0.2 noise
Demo_Gauss.m       % Gaussian noise
Demo_Laplace.m     % Laplace noise
Demo_Sparse.m      % Sparse noise
Demo_Mixture1.m    % Mixture1 noise: Sparse noise + Gaussian noise + Gaussian noise
Demo_Mixture2.m    % Mixture2 noise: EP 0.5 noise + Gaussian noise + Laplace noise

================================================================================================

Main Routine

[label,model,TW,OutU,OutV,llh,llh_BIC,p] = EM_PMoEP(InW,InX,r,param,p,lambda)
%Input:
   InW: d x n x param.k indicator matrices
   InX: d x n input data matrix
   r:   the rank
   param:
      --param.maxiter: maximal iteration number
      --param.OriX: ground truth matrix
      --param.InU: initialized factorized matrice U
      --param.InV: initialized factorized matrice V
      --param.k: the number of mixture components
      --param.display: display the iterative process
      --param.tol: the tolerance for stop
   p: the candidate components
   lambda: the tuning parameter

%Output:
   label: the labels of the noises
   model: model.eta, the precisions of the different EPs
          model.Pi,the mixing coefficients
   W: d x n weighted matrix
   OutU: the final factorized matrix U
   OutV: the final factorized matrix V
   llh:  the log likelihood
   llh_BIC:  the log likelihood used in BIC criterion
   p: the selected components

- =========================================================================    
If you have any quesion, please contact Xiangyong Cao(caoxiangyong45@gmail.com)
