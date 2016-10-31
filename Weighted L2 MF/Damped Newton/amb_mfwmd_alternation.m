function state = amb_mfwmd_alternation(state,m,n,r)

% state = amb_mfwmd_alternation(state[,m,n,r])
%
% Perform one iteration of the alternation algorithm. Provide all
% arguments for initialization of 'state' struct. The fields that are
% required to already exist are: 
%
%   A, B, M, W, L, Lr
%
% Citation: 
%
% See 'amb_mfwmd'

FUNCTION_NAME = 'amb_mfwmd_alternation';

if nargin==4
    % initialize

    state.m = m;
    state.n = n;
    state.I_alt = state.L * eye(r,r);
else
    % minimize
    
	% improve A
	for i=1:state.m
		% the diag matrix can be huge, so must avoid making it
		Wi = state.W(i,:).^2;
		BtWiWi = state.B';
		for j=1:state.n
			BtWiWi(:,j) = BtWiWi(:,j)*Wi(j);  
		end
		state.A(i,:) = ( (BtWiWi*state.B + state.I_alt)\(BtWiWi*state.M(i,:)') )';
	end
    
	% improve B
	for j=1:state.n
		% the diag matrix can be huge, so must avoid making it
		Wj = state.W(:,j).^2;
		AtWjWj = state.A';
		for i=1:state.m
			AtWjWj(:,i) = AtWjWj(:,i)*Wj(i);  
		end
		state.B(j,:) = ( (AtWjWj*state.A + state.I_alt)\(AtWjWj*state.M(:,j)) )';
	end
	
	sumsqrderr = amb_fwmd_errfunc_sumsqrderr(state.A,state.B,state.M,state.W);
	regerr = amb_fwmd_errfunc_regularizeAB(state.A,state.B,state.Lr);
	state.error = sumsqrderr + state.L*regerr;
    
end
    
