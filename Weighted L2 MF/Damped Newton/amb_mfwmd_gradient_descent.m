function state = amb_mfwmd_gradient_descent(state,m,n,r)

% state = amb_mfwmd_gradient_descent(state[,m,n,r])
%
% Perform one iteration of the gradient descent line search algorithm.
% Provide all arguments for initialization of 'state' struct. The fields
% that are required to already exist are: 
%
%   A, B, M, W, L, Lr
%
% Citation: 
%
% See 'amb_mfwmd'

FUNCTION_NAME = 'amb_mfwmd_gradient_descent';

if nargin==4
    % initialize

	state.line_search_tolerance = eps*10;
    state.line_search_limit = 50;

    state.m = m;
    state.n = n;
    state.r = r;
    state.mr = m*r;
    state.nr = n*r;
    
else
    % minimize

	[last_error,fd] = amb_fwmd_errfunc_sumsqrderr(state.A,state.B,state.M,state.W);
	
	[last_regerr,fd_sup] = amb_fwmd_errfunc_regularizeAB(state.A,state.B,state.Lr);
    last_error = last_error + state.L*last_regerr;
	fd = fd + state.L*fd_sup;

	fd = -fd/norm(fd); % want to move downhill

	[x_new,alpha,state.error,success] = amb_fwmd_errfunc_line_search([reshape(state.A',state.mr,1); reshape(state.B',state.nr,1)],fd,state,'sumsqrderr',1,'regularizeAB',state.L,'regratio',state.Lr,'input_checks',0,'line_search_tolerance',state.line_search_tolerance,'line_search_limit',state.line_search_limit);
    
	if ~success
		disp([FUNCTION_NAME ': line search iteration reached its iteration limit.'])  
	end
	
	state.A = reshape(x_new(1:state.mr), state.r, state.m)';
	state.B = reshape(x_new(state.mr+1:end), state.r, state.n)';  
    
end

  
