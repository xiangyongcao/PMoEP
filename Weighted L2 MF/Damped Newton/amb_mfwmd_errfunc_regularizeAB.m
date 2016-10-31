function [error,fd,HA,HB,HC] = amb_mfwmd_errfunc_regularizeAB(A,B,ratio)

%  error = amb_mfwmd_errfunc_regularizeAB(A,B[,ratio])
% [error,first_derivative] = amb_mfwmd_errfunc_regularizeAB(A,B[,ratio])
% [error,first_derivative,hessian] = amb_mfwmd_errfunc_regularizeAB(A,B[,ratio])
% [error,first_derivative,HA,HB,HC] = amb_mfwmd_errfunc_regularizeAB(A,B[,ratio])
%
% error = || A ||^2 + ratio*|| B ||^2 (frobenius norm)
%
% Note: HC is nr-by-r (diagonal blocks 'slid' to the lefthand side)
%       'ratio' is taken as one if omitted.

FUNCTION_NAME = 'amb_mfwmd_errfunc_regularizeAB';

% assume input matrices are consistent
[m,r] = size(A);
n = size(B,1);

if nargin<3
  ratio = 1;
end

if nargout>0
%   disp([FUNCTION_NAME ': calculating error.']) %%% DEBUG

  error = norm(A,'fro')^2 + ratio*norm(B,'fro')^2; % error
  
  if nargout>1
%     disp([FUNCTION_NAME ': calculating first derivative.']) %%% DEBUG

    fd = zeros(m*r+n*r,1);

    % compute first half of d
    for a = 1:m
      for b = 1:r
        index = (a-1)*r+b;   
        fd(index) = 2*A(a,b);
      end
    end
    % compute second half of d
    for c = 1:n
      for d = 1:r
        index = m*r + (c-1)*r+d;
        fd(index) = ratio*2*B(c,d);
      end
    end  
    
    if nargout>2
%       disp([FUNCTION_NAME ': calculating hessian.']) %%% DEBUG

      HA = zeros(m*r,m*r);
      HB = zeros(m*r,n*r);
      HC = zeros(n*r,r);
      
      %    ab cd
      %   /-----\
      % e |  |  |
      % f |  |  |
      %   |--+--| = H
      % g |  |  |
      % h |  |  |
      %   \-----/
      
      % compute upper-left block of H: HA (d2e/dA2)
      for a = 1:m
        for b = 1:r
          second_index = (a-1)*r+b;
          e = a; f = b;
          first_index  = (e-1)*r+f;
          HA(first_index,second_index) = 2;
        end
      end
      
      % compute upper-right block of H: HB (d2e/dAdB)      
      % *** all zeros ***
      
      % compute lower-right block of H: HC (d2e/dB2)
      for c = 1:n
        for d = 1:r
          second_index = d; % (c-1)*r+d;
          g = c; h = d;
          first_index = (g-1)*r+h;
          HC(first_index,second_index) = ratio*2;
        end
      end
      
      if nargout==3
        [d,r] = size(HC);
        HC_full = zeros(d,d);
        for i = 1:r:d
          indices = i:i+r-1;
          HC_full(indices,indices) = HC(indices,:);
        end
        HA = [HA HB; HB' HC_full];
      end
    end
  end
end

% testing shows that the norm(*,'fro')^2 function is faster than
% 1. norm = sum(M(:).^2)
% 2. v=M(:); norm = v'*v;
% 3. norm = trace(M*M');

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CHECK %%%%% CHECK %%%%% CHECK %%%%% CHECK %%%%% CHECK %%%%% CHECK %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout>1
  first_derivative_discrepency = fd - regularizeAB_fd_finite_difference(A,B,ratio); %, pause
  if nargout>2
    [HA_fd,HB_fd,HC_fd] = regularizeAB_H_finite_difference(A,B,ratio);    
    if nargout==3
      hessian_discrepency = HA - [HA_fd HB_fd; HB_fd' HC_fd]; %, pause
    else
      hessian_HA_discrepency = HA - HA_fd; %, pause
      hessian_HB_discrepency = HB - HB_fd; %, pause
      [d,r] = size(HC);
      HC_full = zeros(d,d);
      for i = 1:r:d
        indices = i:i+r-1;
        HC_full(indices,indices) = HC(indices,:);
      end
      hessian_HC_discrepency = HC_full - HC_fd; %, pause
    end
  end
end

if nargout>1
  average_first_derivative_discrepency = sum(first_derivative_discrepency(:))/prod(size(first_derivative_discrepency))
  max_abs_first_derivative_discrepency = max(abs(first_derivative_discrepency(:)))
  if nargout>2
    if nargout==3
      average_hessian_discrepency = sum(hessian_discrepency(:))/prod(size(hessian_discrepency))
      max_abs_hessian_discrepency = max(abs(hessian_discrepency(:)))
    else
      average_hessian_HA_discrepency = sum(hessian_HA_discrepency(:))/prod(size(hessian_HA_discrepency))
      max_abs_hessian_HA_discrepency = max(abs(hessian_HA_discrepency(:)))
      average_hessian_HB_discrepency = sum(hessian_HB_discrepency(:))/prod(size(hessian_HB_discrepency))
      max_abs_hessian_HB_discrepency = max(abs(hessian_HB_discrepency(:)))
      average_hessian_HC_discrepency = sum(hessian_HC_discrepency(:))/prod(size(hessian_HC_discrepency))
      max_abs_hessian_HC_discrepency = max(abs(hessian_HC_discrepency(:)))
    end
  end
end

function fd_fd = regularizeAB_fd_finite_difference(A,B,ratio)
  epsilon = 1e-5;
  
  % assume input matrices are consistent
  [m,r] = size(A);
  n = size(B,1);
      
  fd_fd = zeros(m*r+n*r,1);

  % compute first half of d (de/dA)
  for a = 1:m
    for b = 1:r
      index = (a-1)*r+b;
      
      A_plus = A; A_plus(a,b) = A_plus(a,b) + epsilon;
      plus = amb_mfwmd_errfunc_regularizeAB(A_plus,B,ratio);
      A_minus = A; A_minus(a,b) = A_minus(a,b) - epsilon;
      minus = amb_mfwmd_errfunc_regularizeAB(A_minus,B,ratio);
      
      fd_fd(index) = (plus - minus)/(2*epsilon);
    end
  end
  % compute second half of d (de/DB)
  for c = 1:n
    for d = 1:r
      index = m*r + (c-1)*r+d;
      
      B_plus = B; B_plus(c,d) = B_plus(c,d) + epsilon;
      plus = amb_mfwmd_errfunc_regularizeAB(A,B_plus,ratio);
      B_minus = B; B_minus(c,d) = B_minus(c,d) - epsilon;
      minus = amb_mfwmd_errfunc_regularizeAB(A,B_minus,ratio);
      
      fd_fd(index) = (plus - minus)/(2*epsilon);
    end  
  end
  
return

function [HA_fd,HB_fd,HC_fd] = regularizeAB_H_finite_difference(A,B,ratio)
  epsilon = 1e-5;
  
  % assume input matrices are consistent
  [m,r] = size(A);
  n = size(B,1);
      
  HA_fd = zeros(m*r,m*r);
  HB_fd = zeros(m*r,n*r);
  HC_fd = zeros(n*r,n*r);
  
  %    ab cd
  %   /-----\
  % e |  |  |
  % f |  |  |
  %   |--+--| = H
  % g |  |  |
  % h |  |  |
  %   \-----/
  
  % compute upper-left block of H: HA (d2e/dA2)
  for a = 1:m
    for b = 1:r
      second_index = (a-1)*r+b;
      
      A_plus = A; A_plus(a,b) = A_plus(a,b) + epsilon;
      A_minus = A; A_minus(a,b) = A_minus(a,b) - epsilon;
      
      for e = a:m
        f_start = 1; if e==a, f_start=b; end
        for f = f_start:r
          first_index = (e-1)*r+f;
          
          A_plus_plus = A_plus; A_plus_plus(e,f) = A_plus_plus(e,f) + epsilon; plus_plus = amb_mfwmd_errfunc_regularizeAB(A_plus_plus,B,ratio);
          A_plus_minus = A_plus; A_plus_minus(e,f) = A_plus_minus(e,f) - epsilon; plus_minus = amb_mfwmd_errfunc_regularizeAB(A_plus_minus,B,ratio);
          A_minus_plus = A_minus; A_minus_plus(e,f) = A_minus_plus(e,f) + epsilon; minus_plus = amb_mfwmd_errfunc_regularizeAB(A_minus_plus,B,ratio);
          A_minus_minus = A_minus; A_minus_minus(e,f) = A_minus_minus(e,f) - epsilon; minus_minus = amb_mfwmd_errfunc_regularizeAB(A_minus_minus,B,ratio);
          
          % Abramowitz and Stegun 1972, p. 884
          % (http://gsbwww.uchicago.edu/computing/research/SASManual/ormp/chap5/sect28.htm)
          HA_fd(first_index,second_index) = (plus_plus-plus_minus-minus_plus+minus_minus)/(4*epsilon*epsilon);
        end
      end
    end
  end
  % copy upper into lower
  for i = 1:m*r
    for j = i+1:m*r
      HA_fd(i,j) = HA_fd(j,i);
    end
  end
  
  % compute upper-right block of H: HB (d2e/dAdB)
  for c = 1:n
    for d = 1:r
      second_index = (c-1)*r+d;
      
      B_plus = B; B_plus(c,d) = B_plus(c,d) + epsilon;
      B_minus = B; B_minus(c,d) = B_minus(c,d) - epsilon;
            
      for e = 1:m
        for f = 1:r
          first_index = (e-1)*r+f;

          A_plus = A; A_plus(e,f) = A_plus(e,f) + epsilon;
          A_minus = A; A_minus(e,f) = A_minus(e,f) - epsilon;
          plus_plus = amb_mfwmd_errfunc_regularizeAB(A_plus,B_plus,ratio);
          plus_minus = amb_mfwmd_errfunc_regularizeAB(A_plus,B_minus,ratio);
          minus_plus = amb_mfwmd_errfunc_regularizeAB(A_minus,B_plus,ratio);
          minus_minus = amb_mfwmd_errfunc_regularizeAB(A_minus,B_minus,ratio);
                
          HB_fd(first_index,second_index) = (plus_plus-plus_minus-minus_plus+minus_minus)/(4*epsilon*epsilon);
        end
      end
    end
  end
  
  % compute lower-right block of H: HC (d2e/dB2)
  for c = 1:n
    for d = 1:r
      second_index = (c-1)*r+d;
      
      B_plus = B; B_plus(c,d) = B_plus(c,d) + epsilon;
      B_minus = B; B_minus(c,d) = B_minus(c,d) - epsilon;
      
      for g = c:n
        h_start = 1; if g==c, h_start=d; end
        for h = h_start:r
          first_index = (g-1)*r+h;
          
          B_plus_plus = B_plus; B_plus_plus(g,h) = B_plus_plus(g,h) + epsilon; plus_plus = amb_mfwmd_errfunc_regularizeAB(A,B_plus_plus,ratio);
          B_plus_minus = B_plus; B_plus_minus(g,h) = B_plus_minus(g,h) - epsilon; plus_minus = amb_mfwmd_errfunc_regularizeAB(A,B_plus_minus,ratio);
          B_minus_plus = B_minus; B_minus_plus(g,h) = B_minus_plus(g,h) + epsilon; minus_plus = amb_mfwmd_errfunc_regularizeAB(A,B_minus_plus,ratio);
          B_minus_minus = B_minus; B_minus_minus(g,h) = B_minus_minus(g,h) - epsilon; minus_minus = amb_mfwmd_errfunc_regularizeAB(A,B_minus_minus,ratio);
          
          HC_fd(first_index,second_index) = (plus_plus-plus_minus-minus_plus+minus_minus)/(4*epsilon*epsilon);
        end
      end
    end
  end
  % copy upper into lower
  for i = 1:n*r
    for j = i+1:n*r
      HC_fd(i,j) = HC_fd(j,i);          
    end
  end

return
