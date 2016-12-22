function state = amb_mfwmd_shum(state,m,n,r)

% state = amb_mfwmd_shum(state[,m,n,r])
%
% Perform one iteration of the 'SFM' alternation algorithm. Provide all 
% arguments for initialization of 'state' struct. The fields that are
% required to already exist are: 
%
%   A, B, M, W, L, Lr
%
% Last column of B must be all ones.
%
% Citation: Shum, Ikeuchi and Reddy, "Principal Component Analysis with
% Missing Data ..." PAMI 19(9):855--867, 1995
%
% See 'amb_mfwmd'

FUNCTION_NAME = 'amb_mfwmd_shum';

if nargin==4
    % initialize

	if ~isequal(state.B(:,r), ones(n,1))
		error([FUNCTION_NAME ': Expecting the last row of B to be entirely ones. ''SFM'' solutions only.'])
		%state.B(:,r) = 1;
	end	

    state.m = m;
    state.n = n;
    state.r = r;
    state.r_1 = r-1;
	state.l = ones(n,1); % this variable is called 'el', not 'one'
    state.Ir = state.L * eye(r,r);
	state.I  = state.L * eye(r-1,r-1);
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
		state.A(i,:) = ( (BtWiWi*state.B + state.Ir)\(BtWiWi*state.M(i,:)') )';
	end
    
	% improve B
	P = state.A(:,1:state.r_1);
	t = state.A(:,state.r);
	for j=1:state.n
		% the diag matrix can be huge, so must avoid making it
		Wj = state.W(:,j).^2;
		PtWjWj = P';
		for i=1:state.m
			PtWjWj(:,i) = PtWjWj(:,i)*Wj(i);  
		end
		X(:,j) = (PtWjWj*P + state.I)\(PtWjWj*(state.M(:,j)-t));
	end
	
	state.B = [X' state.l];

	sumsqrderr = amb_fwmd_errfunc_sumsqrderr(state.A,state.B,state.M,state.W);
	regerr = amb_fwmd_errfunc_regularizeAB(state.A,state.B,state.Lr);
	state.error = sumsqrderr + state.L*regerr;
    
  
end

