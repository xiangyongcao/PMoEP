function [m,n,r,valid] = amb_fwmd_matrix_dimension_check(A,B,M,W,FUNCTION_NAME)

% [m,n,r,valid] = amb_fwmd_matrix_dimension_check(A,B,M,W)
%
% Checks consistency of input matrices and returns the key dimensions.
% Correct input is:
%   A is m-by-r
%   B is n-by-r
%   M is m-by-n
%   W is m-by-n

FUNCTION_NAME = 'amb_fwmd_matrix_dimension_check';

valid = 1;

[mA,rA] = size(A);
[nB,rB] = size(B);
[mM,nM] = size(M);
[mW,nW] = size(W);
if mA~=mM, disp([FUNCTION_NAME ': matrix dimensions do not agree - A has a different number of rows to M.']), valid = 0; end
if nB~=nM, disp([FUNCTION_NAME ': matrix dimensions do not agree - B'' has a different number of columns to M.']), valid = 0; end
if rA~=rB, disp([FUNCTION_NAME ': matrix dimensions do not agree - A has a different number of columns to B.']), valid = 0; end
if mW~=mM, disp([FUNCTION_NAME ': matrix dimensions do not agree - W has a different number of rows to M.']), valid = 0; end
if nW~=nM, disp([FUNCTION_NAME ': matrix dimensions do not agree - W has a different number of columns to M.']), valid = 0; end

m = mM;
n = nM;
r = rA;

