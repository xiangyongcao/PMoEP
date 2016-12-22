function [x1,x2] = amb_block_matrix_solve(A,B,C,b1,b2,r,C_not_block_diagonal)

% [x1,x2] = amb_hessian_solve_sparse(A11,A12,A22,b1,b2,r[,A22_not_block_diagonal])
%
% Solves the matrix equation Ax = b when A is symmetric and of the form
% where A22 has non-zero values only in diagonal r-by-r blocks (unless the
% last optional argument is set to 1).
% 
%                    / A11 | A12 \ / x1 \   / b1 \
%                    |-----+-----| |----| = |----|
%                    \ A12'| A22 / \ x2 /   \ b2 /
%
% Efficiency wise, it is assumed that the dimension of A11 is smaller than
% that of A22. Also, A22 can be given as a nr-by-r rectangular matrix (i.e.
% all the r-by-r diagonal blocks slid to the left hand side of the matrix).

FUNCTION_NAME = 'amb_block_matrix_solve';

% TO DO: A11 and A22 are symetrical and so only half of either needs to be
% provided.

% In an attempt to avoid confusion, more varied names are used for the
% blocks of the matrix to be inverted:
%
%                    /  A  |  B  \ / x1 \   / b1 \
%                    |-----+-----| |----| = |----|
%                    \  B' |  C  / \ x2 /   \ b2 /
%
% and the inverse equation:
%
%                    / x1 \   /  P  |  Q  \ / b1 \
%                    |----| = |-----+-----| |----|
%                    \ x2 /   \  Q' |  R  / \ b2 /
%

if nargin<7
  C_not_block_diagonal = 0;
end

nr = size(C,1);
mr = size(B,1);
n = floor(nr/r);

if nr<mr, disp([FUNCTION_NAME ': warning: A11 is larger than A22.']), end

C_square = 1;
if size(C,2)==r, C_square = 0; end

CinvBt = zeros(nr,mr);
Cinvb2 = zeros(nr,1);
if C_not_block_diagonal
  CinvBt = C\B';
  Cinvb2 = C\b2;
else
  for i = 1:n
    offset = (i-1)*r;
    if C_square==1
      Ci = C(offset+1:offset+r,offset+1:offset+r); 
    else
      Ci = C(offset+1:offset+r,1:r);
    end
    
    Cii = inv(Ci);
    tmp = Cii*B(:,offset+1:offset+r)';
    CinvBt(offset+1:offset+r,:) = tmp;
    Cinvb2(offset+1:offset+r) = Cii*b2(offset+1:offset+r);
    % the following possibly cheaper for large matrices...?
    %tmp = Ci\[B(:,offset+1:offset+r)' b2(offset+1:offset+r)];
    %CinvBt(offset+1:offset+r,:) = tmp(:,1:end-1);
    %Cinvb2(offset+1:offset+r) = tmp(:, end);
  end
end

% taking ^ to mean inverse and noting that C^' = C^
% P = (A - BC^B')^
% Q = -P'BC^ = -(A - BC^B')^'BC^
% R = C^(I - B'Q) = C^(I + B'(A - BC^B')^'BC^)

spCinvBt = sparse(CinvBt);
% the following should be set up outside this function
spA = sparse(A);
spB = sparse(B);

% assume m<n, therefore P smaller than R
Pinv = spA - spB*spCinvBt; % this line is twice as quick with sparse
%Pinv_b1 = Pinv\b1;
%Pinvt_B_Cinv_b2 = Pinv\(B*Cinvb2); % Pinv = Pinv'
%x1 =   Pinv_b1 - Pinvt_B_Cinv_b2; 

x1 = Pinv\(b1 - B*Cinvb2);
x2 = Cinvb2 -  CinvBt*x1; 

x1 = full(x1);
x2 = full(x2);