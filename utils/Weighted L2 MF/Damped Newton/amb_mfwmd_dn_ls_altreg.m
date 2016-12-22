function state = amb_mfwmd_dn_line_search(state,m,n,r)

% state = amb_mfwmd_alternation(state[,m,n,r])
%
% Perform one iteration of the damped Newton line search algorithm. Provide
% all arguments for initialization of 'state' struct. The fields that are 
% required to already exist are: 
%
%   A, B, M, W, L, Lr
%
% Citation: 
%
% See 'amb_mfwmd'

FUNCTION_NAME = 'amb_mfwmd_dn_line_search';

if nargin==4
    % initialize

	state.line_search_tolerance = 1e-4;
    state.line_search_limit = 20;

    state.m = m;
    state.n = n;
    state.r = r;
    state.mr = m*r;
    state.nr = n*r;
    state.mn = m*n;
    
	state.I_dn = eye(m*r,m*r);
else
    % minimize

	[last_error,fd,HA,HB,HC] = amb_fwmd_errfunc_sumsqrderr(state.A,state.B,state.M,state.W);
	
	[last_regerr,fd_sup,HA_sup,HB_sup,HC_sup] = amb_fwmd_errfunc_regularizeAB(state.A,state.B,state.Lr);
    last_error = last_error + state.L*last_regerr;
	fd = fd + state.L*fd_sup;
	HA = HA + state.L*HA_sup;
	HB = HB + state.L*HB_sup;
	HC = HC + state.L*HC_sup;
 
	% "altreg"
	% make it resort to a more coordinate-descent-like direction rather
	% than a gradient descent direction
	HA = HA*(state.lambda+1);
	HC = HC*(state.lambda+1);
	
    % to ensure invertability, regularize as usual
    HESSIAN_STABILITY = 1e-7;
	HA = HA + HESSIAN_STABILITY*state.I_dn;
	% HCI = HC + lambda*eye(nr);
	% above line implemented as [
	for i=1:state.r
		indices = [i:state.r:state.nr] + (i-1)*state.nr;
		HC(indices) = HC(indices) + HESSIAN_STABILITY;
	end
	% ] Note: C is not square!
    
	[h_mr,h_nr] = amb_block_matrix_solve(HA,HB,HC,fd(1:state.mr),fd(state.mr+1:end),state.r);
	h = -[h_mr;h_nr];

	[x_new,alpha,state.error,success] = amb_fwmd_errfunc_line_search([reshape(state.A',state.mr,1); reshape(state.B',state.nr,1)],h,state,'sumsqrderr',1,'regularizeAB',state.L,'regratio',state.Lr,'input_checks',0,'line_search_tolerance',state.line_search_tolerance,'line_search_limit',state.line_search_limit);
  
	if ~success
		disp([FUNCTION_NAME ': line search iteration reached its iteration limit (probably)'])  
	end
	
	%%% update lambda using alpha
	% lambda = lambda/sqrt(alpha); % this seems to react too slowly and end up with almost linear converge
	state.lambda = state.lambda/alpha;

	state.A = reshape(x_new(1:state.mr), state.r, state.m)';
	state.B = reshape(x_new(state.mr+1:end), state.r, state.n)';  
end

