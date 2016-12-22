function state = amb_mfwmd_alt_dn_hybrid_v02(state,m,n,r,record_switches)

% state = amb_mfwmd_alt_dn_hybrid_v02(state[,m,n,r,record_switches])
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
% *Another combination of alternation and 'Damped Newton'. This algorithm
% performs a set number of alternation steps and then goes into 'hybrid'
% mode. After any iteration, if damped Newton's lambda is too large then a
% differnent number of alternation steps are performed. 
%
% See 'amb_mfwmd'

FUNCTION_NAME = 'amb_mfwmd_alt_dn_hybrid_v02';

if nargin==5
    % initialize

    state.record_switches = record_switches;
    if record_switches, state.switches = [0]; end
    
    state.alternation_countdown = 8; % the number of alternation steps before starting damped Newton
    state.alternation_countdown_reset = 8; % the number of alternation steps when lambda gets bad

    state.lambda_max = 1e3;
    
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
        state = amb_mfwmd_damped_newton(state);
        if state.record_switches, state.switches(end+1) = 0; end
        if state.lambda > state.lambda_max
            state.lambda = state.lambda_max;
            state.alternation_countdown = state.alternation_countdown_reset; 
        end
        if state.alternation_countdown<0, state.alternation_countdown = state.alternation_countdown + 1; end
    end
  
end


