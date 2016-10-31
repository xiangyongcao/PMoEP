function state = amb_mfwmd_alt_dn_hybrid_v03(state,m,n,r,record_switches)

% state = amb_mfwmd_alt_dn_hybrid_v03(state[,m,n,r,record_switches])
%
% Perform one iteration of the alternation/damped Newton hybrid (version 3)
% algorithm*. Provide all arguments for initialization of 'state' struct.
% The fields that are required to already exist are:  
%
%   A, B, M, W, L, Lr, lambda, lambda_step
%
% Citation: Buchanan and Fitzgibbon, "Damped Newton Algorithms for Matrix
% Factorization with Missing Data" CVPR'05
%
% *Another combination of alternation and damped Newton. This algorithm
% performs a set number of alternation steps and then goes into 'hybrid'
% mode. In each iteration the Levenberg-Marquardt style loop is entered,
% but the hessian-regularizing lambda (that controls step length and
% gradient-descent similarity) is monitored. If it becomes too large or is
% worse than the last iteration then a differnent number of alternation
% steps are performed.
%
% See 'amb_mfwmd'

FUNCTION_NAME = 'amb_mfwmd_alt_dn_hybrid_v03';

if nargin==5
    % initialize

    state.record_switches = record_switches;
    if state.record_switches, state.switches = [0]; end
    
    state.alternation_countdown = 8; % the number of alternation steps before starting damped Newton
    state.alternation_countdown_reset = 0; % the number of EXTRA alternation steps when lambda gets bad

    state.lambda_max = 1e4;
    
    state.m = m;
    state.n = n;
    state.r = r;

    state.mr = m*r;
    state.nr = n*r;
    state.mn = m*n;
    
	state.I_dn = eye(m*r,m*r);
	state.I_alt = state.L*eye(r,r);
	
else
    % minimize

    if state.alternation_countdown>0
    	%%%% Alternation
        state = amb_mfwmd_alternation(state);
        state.alternation_countdown = state.alternation_countdown - 1;
        if state.record_switches, state.switches(end+1) = 1; end
    else
    	%%%% Damped Newton
        if state.record_switches, state.switches(end+1) = 0; end % will be overridden if alternation used
        
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
		
        last_lambda = state.lambda;
		state.lambda = state.lambda/(state.lambda_step^2);
		
		improvement = -inf;
		while improvement<-state.tolerance
			state.lambda = state.lambda*state.lambda_step;
			if state.lambda>1e12, disp([FUNCTION_NAME ': step size too small; giving up.']), break, end
			
            if state.lambda > state.lambda_step*last_lambda || state.lambda > state.lambda_max; % if above initial level or maxed out
                state = amb_mfwmd_alternation(state);
                state.alternation_countdown = state.alternation_countdown_reset;
                if state.record_switches, state.switches(end) = 1; end
				state.lambda = state.lambda*state.lambda_step; % keep lambda progression consistent
                break
            else
				HAI = HA + state.lambda*state.I_dn;
				% HCI = HC + lambda*eye(nr);
				% above line implemented as [
				HCI = HC;
				for i=1:state.r
					indices = [i:state.r:state.nr] + (i-1)*state.nr;
					HCI(indices) = HCI(indices)+ state.lambda;
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
		end % while
        
        if state.alternation_countdown<0, state.alternation_countdown = state.alternation_countdown + 1; end
    end
  
end


