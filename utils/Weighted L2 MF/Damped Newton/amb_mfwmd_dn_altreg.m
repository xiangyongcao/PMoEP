function state = amb_mfwmd_damped_newton(state,m,n,r)

% state = amb_mfwmd_damped_newton(state[,m,n,r])
%
% Perform one iteration of the damped Newton algorithm with 'alternation
% style' stabalization for the hessian. Provide all arguments for
% initialization of 'state' struct. The fields that are required to already
% exist are:  
%
%   A, B, M, W, tolerance, lambda, lambda_step 
%
% Citation: Buchanan and Fitzgibbon, "Damped Newton Algorithms for Matrix
% Factorization with Missing Data" CVPR'05
%
% See 'amb_mfwmd'

FUNCTION_NAME = 'amb_mfwmd_damped_newton';

if nargin==4
    % initialize

    state.m = m;
    state.n = n;
    state.r = r;
    state.mr = m*r;
    state.nr = n*r;
    state.mn = m*n;

    state.I_dn = eye(state.mr);
else
    % minimize
    
	[last_error,fd,HA,HB,HC] = amb_fwmd_errfunc_sumsqrderr(state.A,state.B,state.M,state.W);
	
	[last_regerr,fd_sup,HA_sup,HB_sup,HC_sup] = amb_fwmd_errfunc_regularizeAB(state.A,state.B,state.Lr);
    last_error = last_error + state.L*last_regerr;
	fd = fd + state.L*fd_sup;
	HA = HA + state.L*HA_sup;
	HB = HB + state.L*HB_sup;
	HC = HC + state.L*HC_sup;
    
	vec_A = reshape(state.A',state.mr,1);
	vec_B = reshape(state.B',state.nr,1);
	
	fd_mr = fd(1:state.mr);
	fd_nr = fd(state.mr+1:end);
	
	state.lambda = state.lambda/(state.lambda_step^2);
	
	improvement = -inf;
	while improvement<-state.tolerance
		state.lambda = state.lambda*state.lambda_step;
		if state.lambda>1e12, disp([FUNCTION_NAME ': step size too small; giving up.']), break, end
		
		% "altreg"
		% make it resort to a more coordinate-descent-like direction rather
		% than a gradient descent direction
		HA = HA*(state.lambda+1);
		HC = HC*(state.lambda+1);
		
		% to ensure invertability, regularize as usual
        HESSIAN_STABILITY = 1e-9;
		HAI = HA + HESSIAN_STABILITY*state.I_dn;
		% HCI = HC + lambda*eye(nr);
		% above line implemented as [
		HCI = HC;
		for i=1:state.r
			indices = [i:state.r:state.nr] + (i-1)*state.nr;
			HCI(indices) = HCI(indices)+ HESSIAN_STABILITY;
		end
		% ] Note: C is not square!
        
		[h_mr,h_nr] = amb_block_matrix_solve(HAI,HB,HCI,fd_mr,fd_nr,state.r);
		
		vec_A_new = vec_A - h_mr; 
		vec_B_new = vec_B - h_nr;
		
		state.A = reshape(vec_A_new,state.r,state.m)';
		state.B = reshape(vec_B_new,state.r,state.n)';
		state.error = amb_fwmd_errfunc_sumsqrderr(state.A,state.B,state.M,state.W);
		state.error = state.error + state.L * amb_fwmd_errfunc_regularizeAB(state.A,state.B,state.Lr);
		improvement = last_error - state.error; % expecting positive improvement      
	end
  
end