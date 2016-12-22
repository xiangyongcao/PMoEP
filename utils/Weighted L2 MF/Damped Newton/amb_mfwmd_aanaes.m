function state = amb_mfwmd_aanaes(state,m,n,r)

% state = amb_mfwmd_aanaes(state[,m,n,r])
%
% Perform one iteration of the Aanaes et al alternation algorithm. Provide
% all arguments for initialization of 'state' struct. The fields that are
% required to already exist are: 
%
%   A, B, M, W, L, Lr
%
% Citation: Aanaes, Fisker, Astrom and Carstensen, "Robust Factorization",
% PAMI, 24(9):1215-1225, 2002.
%
% See 'amb_mfwmd'

FUNCTION_NAME = 'amb_mfwmd_aanaes';

IMPUTATION_INITIALIZATION = 1;

if nargin==4
    % initialize

    state.m = m;
    state.n = n;
    state.r = r;
    state.I = state.L * eye(r,r);
    
    %%%%% set up initial guess for M
	if IMPUTATION_INITIALIZATION
		state.Mhat = state.M;
		Mx = state.M(1:2:end,:); Wx = state.W(1:2:end,:);
		My = state.M(2:2:end,:); Wy = state.W(2:2:end,:);
		if mod(m,2)==0 && isequal(Wx,Wy)
			% find average visible coords (assumes tracking measurement matrix)
			disp([FUNCTION_NAME ': assuming vertically alternating x-y coordinates in M for fill averages.'])
			MWx = Wx.*Mx; MWy = Wy.*My;
			x_avg = sum(MWx(:))/sum(Wx(:));
			y_avg = sum(MWy(:))/sum(Wy(:));
			% fill unknown coord with weighted average
			Mx(Wx==0) = x_avg;
			My(Wy==0) = y_avg;
			state.Mhat(1:2:end) = Mx;
			state.Mhat(2:2:end) = My;
			clear MWx MWy
		else 
			% find average visible element
			disp([FUNCTION_NAME ': assuming general matrix for fill average.'])
			MW = state.W.*state.M;
			average = sum(MW(:))/sum(state.W(:));
			% fill unknown elements with weighted average
			state.Mhat(state.W==0) = average;
			clear MW
		end
		clear Mx My Wx Wy
	else
		state.Mhat = amb_merge_matrices(state.W,state.M,state.A*state.B');  
	end
	

else
    % minimize
     
	% estimate A
	[U,S,V] = svd(state.Mhat,0);
	state.A = U(:,1:state.r);
	
	% improve B
	for j=1:state.n
		% the diag matrix can be huge, so must avoid making it
		Wj = state.W(:,j).^2;
		AtWjWj = state.A';
		for i=1:state.m
			AtWjWj(:,i) = AtWjWj(:,i)*Wj(i);  
		end
		state.B(j,:) = ( (AtWjWj*state.A + state.I)\(AtWjWj*state.Mhat(:,j)) )';
	end
	
	state.Mhat = amb_merge_matrices(state.W,state.M,state.A*state.B');

	sumsqrderr = amb_fwmd_errfunc_sumsqrderr(state.A,state.B,state.M,state.W);
	regerr = amb_fwmd_errfunc_regularizeAB(state.A,state.B,state.Lr);
	state.error = sumsqrderr + state.L*regerr;
    
end