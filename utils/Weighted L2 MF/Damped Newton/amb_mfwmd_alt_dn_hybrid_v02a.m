function state = amb_mfwmd_alt_dn_hybrid_v02a(state,m,n,r,record_switches)

% state = amb_mfwmd_alt_dn_hybrid_v02a(state[,m,n,r,record_switches])
%
% Perform one iteration of the alternation/damped Newton hybrid (version 2a)
% algorithm*. Provide all arguments for initialization of 'state' struct.
% The fields that are required to already exist are:  
%
%   A, B, M, W, L, Lr, lambda
%
% Citation: Buchanan and Fitzgibbon, "Damped Newton Algorithms for Matrix
% Factorization with Missing Data" CVPR'05
%
% *Another combination of alternation and damped Newton with line search.
% This algorithm performs a set number of alternation steps and then goes
% into 'hybrid' mode. After any iteration, if damped Newton's lambda is too
% large then a differnent number of alternation steps are performed. 
%
% See 'amb_mfwmd'

FUNCTION_NAME = 'amb_mfwmd_alt_dn_hybrid_v02a';

if nargin==5
    % initialize

    state.record_switches = record_switches;
    if record_switches, state.switches = [0]; end
    
    state.alternation_countdown = 8; % the number of alternation steps before starting damped Newton
    state.alternation_countdown_reset = 8; % the number of alternation steps when lambda gets bad

    state.lambda_max = 1e3;
    
    state.line_search_tolerance = 1e-4;
    state.line_search_limit = 20;

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
    	%%%% Damped Newton Line Search
        state = amb_mfwmd_dn_line_search(state);
        if state.record_switches, state.switches(end+1) = 0; end
        if state.lambda > state.lambda_max
            state.lambda = state.lambda_max;
            state.alternation_countdown = state.alternation_countdown_reset; 
        end
        if state.alternation_countdown<0, state.alternation_countdown = state.alternation_countdown + 1; end
    end
  
end


