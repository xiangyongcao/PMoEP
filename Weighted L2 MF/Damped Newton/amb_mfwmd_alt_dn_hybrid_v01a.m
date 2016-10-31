function state = amb_mfwmd_alt_dn_hybrid_v01a(state,m,n,r,record_switches)

% state = amb_mfwmd_alt_dn_hybrid_v01a(state[,m,n,r,record_switches])
%
% Perform one iteration of the alternation/damped Newton hybrid (version 1a)
% algorithm*. Provide all arguments for initialization of 'state' struct.
% The fields that are required to already exist are:  
%
%   A, B, M, W, L, Lr, lambda
%
% Citation: Buchanan and Fitzgibbon, "Damped Newton Algorithms for Matrix
% Factorization with Missing Data" CVPR'05
%
% *A combination of alternation and damped Newton with line search: in each
% iteration, calculates new points in error-space for both alternation and
% damped Newton line search and picks the best of the two.
%
% See 'amb_mfwmd'

FUNCTION_NAME = 'amb_mfwmd_alt_dn_hybrid_v01a';

if nargin==5
    % initialize

    state.record_switches = record_switches;
    if record_switches, state.switches = [0]; end
    
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

	%%%% Alternation    
    alt_state = amb_mfwmd_alternation(state);
      
	%%%% Damped Newton Line Search
    dnls_state = amb_mfwmd_dn_line_search(state);
    
    if alt_state.error < dnls_state.error
        state = alt_state;
        state.lambda = dnls_state.lambda;
        if state.record_switches, state.switches(end+1) = 1; end
    else
        state = dnls_state;
        if state.record_switches, state.switches(end+1) = 0; end
    end
  
end

