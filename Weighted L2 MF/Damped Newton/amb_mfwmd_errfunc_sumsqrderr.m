function [error,fd,HA,HB,HC] = amb_mfwmd_errfunc_sumsqrderr(A,B,M,W)

%  error = amb_mfwmd_errfunc_sumsqrderr(A,B,M,W)
% [error,first_derivative] = amb_mfwmd_errfunc_sumsqrderr(A,B,M,W)
% [error,first_derivative,hessian] = amb_mfwmd_errfunc_sumsqrderr(A,B,M,W)
% [error,first_derivative,HA,HB,HC] = amb_mfwmd_errfunc_sumsqrderr(A,B,M,W)
%
% error = || W.*( M - A*B' ) ||^2 (frobenius norm)
%
% Note: HC is nr-by-r (diagonal blocks 'slid' to the lefthand side)

FUNCTION_NAME = 'amb_mfwmd_errfunc_sumsqrderr';

% assume input matrices are consistent
[m,r] = size(A);
n = size(B,1);
      
USE_MEX = 1;
if USE_MEX
    % the mex routines cannot deal with sparse matrices
    if issparse(A), A = full(A); end
    if issparse(B), B = full(B); end
    if issparse(M), M = full(M); end
    if issparse(W), W = full(W); end
end

if nargout>0
%   disp([FUNCTION_NAME ': calculating error.']) %%% DEBUG

  error = norm(W.*(M-A*B'),'fro')^2; % error
  
  if nargout>1
%     disp([FUNCTION_NAME ': calculating first derivative.']) %%% DEBUG

    if USE_MEX
      fd = amb_mfwmd_errfunc_sumsqrderr_gradient_mex(A,B,M,W); % first derivative
    else
      fd = zeros(m*r+n*r,1);
      
      % compute first half of d
      for a = 1:m
        for b = 1:r
          index = (a-1)*r+b;   
          for j = 1:n
            AapBjp_sum = 0;
            for p=1:r
              AapBjp_sum = AapBjp_sum + A(a,p)*B(j,p);  
            end
            fd(index) = fd(index) + W(a,j)*W(a,j)*B(j,b)*(AapBjp_sum - M(a,j));
          end
          fd(index) = fd(index)*2;
        end
      end
      % compute second half of d
      for c = 1:n
        for d = 1:r
          index = m*r + (c-1)*r+d;
          for i = 1:m
            AipBcp_sum = 0;
            for p=1:r
              AipBcp_sum = AipBcp_sum + A(i,p)*B(c,p);  
            end
            fd(index) = fd(index) + W(i,c)*W(i,c)*A(i,d)*(AipBcp_sum - M(i,c));
          end
          fd(index) = fd(index) * 2;
        end
      end  
    end
    
    if nargout>2
%       disp([FUNCTION_NAME ': calculating hessian.']) %%% DEBUG

      if USE_MEX
        [HA,HB,HC] = amb_mfwmd_errfunc_sumsqrderr_hessian_mex(A,B,M,W); % hessian
      else
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
        
        % compute upper-left block of H: HA
        for a = 1:m
          for b = 1:r
            second_index = (a-1)*r+b;
            for f = 1:r
              first_index = (a-1)*r+f;
              for j = 1:n
                HA(first_index,second_index) = HA(first_index,second_index) + W(a,j)*W(a,j)*B(j,b)*B(j,f);
              end
              HA(first_index,second_index) = HA(first_index,second_index) * 2;
            end
          end
        end
        
        % compute upper-right block of H: HB
        for c = 1:n
          for d = 1:r
            second_index = (c-1)*r+d;
            for e = 1:m
              for f = 1:r
                first_index = (e-1)*r+f;
                AepBcp_sum = 0;
                if d==f
                  for p = 1:r
                    AepBcp_sum = AepBcp_sum + A(e,p)*B(c,p);  
                  end
                  AepBcp_sum = AepBcp_sum - M(e,c);  
                end
                HB(first_index,second_index) = 2*W(e,c)*W(e,c)*(A(e,d)*B(c,f) + AepBcp_sum); 
                % H(first_index,second_index) = 2*W(e,c)*W(e,c)*(A(e,f)*B(c,d) + AepBcp_sum); % transpose: d<>f (wrong I reckon)
              end
            end
          end
        end
        
        % compute lower-right block of H: HC
        for c = 1:n
          for d = 1:r
            second_index = d; % (c-1)*r+d;
            for h = 1:r
              first_index = (c-1)*r+h;
              for i = 1:m
                HC(first_index,second_index) = HC(first_index,second_index) + W(i,c)*W(i,c)*A(i,d)*A(i,h);
              end
              HC(first_index,second_index) = HC(first_index,second_index) * 2;
            end
          end
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
  first_derivative_discrepency = fd - sumsqrderr_fd_finite_difference(A,B,M,W); %, pause
  if nargout>2
    [HA_fd,HB_fd,HC_fd] = sumsqrderr_H_finite_difference(A,B,M,W);    
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

function fd_fd = sumsqrderr_fd_finite_difference(A,B,M,W)
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
      plus = amb_mfwmd_errfunc_sumsqrderr(A_plus,B,M,W);
      A_minus = A; A_minus(a,b) = A_minus(a,b) - epsilon;
      minus = amb_mfwmd_errfunc_sumsqrderr(A_minus,B,M,W);
      
      fd_fd(index) = (plus - minus)/(2*epsilon);
    end
  end
  % compute second half of d (de/DB)
  for c = 1:n
    for d = 1:r
      index = m*r + (c-1)*r+d;
      
      B_plus = B; B_plus(c,d) = B_plus(c,d) + epsilon;
      plus = amb_mfwmd_errfunc_sumsqrderr(A,B_plus,M,W);
      B_minus = B; B_minus(c,d) = B_minus(c,d) - epsilon;
      minus = amb_mfwmd_errfunc_sumsqrderr(A,B_minus,M,W);
      
      fd_fd(index) = (plus - minus)/(2*epsilon);
    end  
  end
  
return

function [HA_fd,HB_fd,HC_fd] = sumsqrderr_H_finite_difference(A,B,M,W)
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
          
          A_plus_plus = A_plus; A_plus_plus(e,f) = A_plus_plus(e,f) + epsilon; plus_plus = amb_mfwmd_errfunc_sumsqrderr(A_plus_plus,B,M,W);
          A_plus_minus = A_plus; A_plus_minus(e,f) = A_plus_minus(e,f) - epsilon; plus_minus = amb_mfwmd_errfunc_sumsqrderr(A_plus_minus,B,M,W);
          A_minus_plus = A_minus; A_minus_plus(e,f) = A_minus_plus(e,f) + epsilon; minus_plus = amb_mfwmd_errfunc_sumsqrderr(A_minus_plus,B,M,W);
          A_minus_minus = A_minus; A_minus_minus(e,f) = A_minus_minus(e,f) - epsilon; minus_minus = amb_mfwmd_errfunc_sumsqrderr(A_minus_minus,B,M,W);
          
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
          plus_plus = amb_mfwmd_errfunc_sumsqrderr(A_plus,B_plus,M,W);
          plus_minus = amb_mfwmd_errfunc_sumsqrderr(A_plus,B_minus,M,W);
          minus_plus = amb_mfwmd_errfunc_sumsqrderr(A_minus,B_plus,M,W);
          minus_minus = amb_mfwmd_errfunc_sumsqrderr(A_minus,B_minus,M,W);
                
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
          
          B_plus_plus = B_plus; B_plus_plus(g,h) = B_plus_plus(g,h) + epsilon; plus_plus = amb_mfwmd_errfunc_sumsqrderr(A,B_plus_plus,M,W);
          B_plus_minus = B_plus; B_plus_minus(g,h) = B_plus_minus(g,h) - epsilon; plus_minus = amb_mfwmd_errfunc_sumsqrderr(A,B_plus_minus,M,W);
          B_minus_plus = B_minus; B_minus_plus(g,h) = B_minus_plus(g,h) + epsilon; minus_plus = amb_mfwmd_errfunc_sumsqrderr(A,B_minus_plus,M,W);
          B_minus_minus = B_minus; B_minus_minus(g,h) = B_minus_minus(g,h) - epsilon; minus_minus = amb_mfwmd_errfunc_sumsqrderr(A,B_minus_minus,M,W);
          
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
