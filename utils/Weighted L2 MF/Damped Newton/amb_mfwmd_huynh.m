function state = amb_mfwmd_huynh(state,m,n,r)

% state = amb_mfwmd_huynh(state[,m,n,r])
%
% Perform one iteration of the Huynh et al. alternation algorithm. Provide
% all arguments for initialization of 'state' struct. The fields that are 
% required to already exist are: 
%
%   A, B, M, W
%
% Last column of B must be all ones.
%
% Citation: Huynh, Hartley and Heyden, "Outlier correction in image
% sequences for the affine camera" ICCV'03
%
% See 'amb_mfwmd'

FUNCTION_NAME = 'amb_mfwmd_huynh';

IMPUTATION_INITIALIZATION = 1;

if nargin==4
    % initialize

	if ~isequal(state.B(:,r), ones(n,1))
		error([FUNCTION_NAME ': Expecting the last row of B to be entirely ones. ''SFM'' solutions only.'])
		%state.B(:,r) = 1;
	end	

	%%%%% set up initial guess for M
	if IMPUTATION_INITIALIZATION
		state.Mtilde = state.M;
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
			state.Mtilde(1:2:end) = Mx;
			state.Mtilde(2:2:end) = My;
			clear MWx MWy
		else 
			% find average visible element
			disp([FUNCTION_NAME ': assuming general matrix for fill average.'])
			MW = state.W.*state.M;
			average = sum(MW(:))/sum(state.W(:));
			% fill unknown elements with weighted average
			state.Mtilde(state.W==0) = average;
			clear MW
		end
		clear Mx My Wx Wy
	else
		state.Mtilde = amb_merge_matrices(state.W,state.M,state.A*state.B');  
	end
	
	state.l = ones(n,1); % this variable is called 'el', not 'one'
	state.X = nan*zeros(r-1,n);
    
	state.c = 1.5;             % correction constant
	state.vs = m*ones(n,1);    % v stores the number of unmodified elements of the jth column of M
	state.sigmas = zeros(n,1); % observation variances
	state.S = zeros(m,n);      % matrix of 'standard errors'
	
	state.m = m;
    state.n = n;
    state.r = r;
    
    warning([FUNCTION_NAME ': the Huynh alternation scheme is not suitable for missing data problems.']);
    warning([FUNCTION_NAME ': the Huynh alternation scheme does not minimize the regularized error function; only sum-squared error is returned.']);
else
    % minimize

	% improve A
	state.A = ((state.B'*state.B)\state.B'*state.Mtilde')'; % equiv  = Mtilde*B*inv(B'*B);
	P = state.A(:,1:state.r-1);
	t = state.A(:,state.r);
	
	% calculate residuals
	H = P*((P'*P)\P');
	Mhat = H*(state.Mtilde - t*state.l') + t*state.l';
	R = state.Mtilde - Mhat; % residuals
	
	% update sigmas
	for j=1:state.n
        if state.vs(j)==0, state.sigmas(j) = 1e9; % large number
        else state.sigmas(j) = sum(R(:,j))/( (state.m-(state.r-1))*(state.vs(j)/state.m)^2  ); % observation variances
        end
	end
	
	% calculate 'standard errors'
	hs = diag(H);
	for i=1:state.m
		sqrt_1_hi = sqrt(1-hs(i));
		for j=1:state.n
			state.S(i,j) = sqrt_1_hi*state.sigmas(j);
		end
	end
	
	% update M^\tilde (combine code in above loops)
	state.vs = state.m*ones(state.n,1);
	state.S = state.c*state.S;
	for i=1:state.m
		for j=1:state.n
			if R(i,j) < - state.S(i,j)
				state.Mtilde(i,j) = Mhat(i,j) - state.S(i,j);
				state.vs(j) = state.vs(j) - 1;
			elseif R(i,j) > state.S(i,j)
				state.Mtilde(i,j) = Mhat(i,j) + state.S(i,j);
				state.vs(j) = state.vs(j) - 1;
			end
		end
	end
	
	% improve B
	state.X = (P'*P)\(P'*(state.Mtilde-t*state.l')); % equiv  = ( Mtilde' - l*t')*P*inv(P'*P) )';
	state.B = [state.X' state.l];
	
	state.Mtilde = state.A*state.B';
	
	state.error = amb_fwmd_errfunc_sumsqrderr(state.A,state.B,state.M,state.W);
  
end

