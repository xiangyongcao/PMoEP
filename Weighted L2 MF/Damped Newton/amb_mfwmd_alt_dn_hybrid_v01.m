function state = amb_mfwmd_alt_dn_hybrid_v01(state,m,n,r,record_switches)

% state = amb_mfwmd_alt_dn_hybrid_v01(state[,m,n,r,record_switches])
%
% Perform one iteration of the alternation/damped Newton hybrid (version 1)
% algorithm*. Provide all arguments for initialization of 'state' struct.
% The fields that are required to already exist are:  
%
%   A, B, M, W, L, Lr, lambda, lambda_step
%
% Citation: Buchanan and Fitzgibbon, "Damped Newton Algorithms for Matrix
% Factorization with Missing Data" CVPR'05
%
% *A combination of alternation and damped Newton: in each iteration,
% calculates new points in error-space for both alternation and damped
% Newton and picks the best of the two.
%
% See 'amb_mfwmd'

FUNCTION_NAME = 'amb_mfwmd_alt_dn_hybrid_v01';

if nargin==5
    % initialize

    state.record_switches = record_switches;
    if record_switches, state.switches = [0]; end
    
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
      
	%%%% Damped Newton 
    dn_state = amb_mfwmd_damped_newton(state);
    
    if alt_state.error < dn_state.error
        state = alt_state;
        state.lambda = dn_state.lambda;
        if state.record_switches, state.switches(end+1) = 1; end
    else
        state = dn_state;
        if state.record_switches, state.switches(end+1) = 0; end
    end
  
end

