function state = amb_mfwmd_brandt(state,m,n,r)

% state = amb_mfwmd_brandt(state[,m,n,r])
%
% Perform one iteration of the Brandt alternation algorithm. Provide all
% arguments for initialization of 'state' struct. The fields that are
% required to already exist are: 
%
%   A, B, M, W
%
% Last column of B must be all ones.
%
% Citation: Brandt, "Closed-form solutions for affine reconstruction under
% missing data" ECCV'02 Workshop.
%
% See 'amb_mfwmd'

FUNCTION_NAME = 'amb_mfwmd_brandt';

if nargin==4
    % initialize

	% TO DO: include a 'transpose' option so that the following warning may be
	% heeded easily. At the moment, the code relies on the tracks as columns,
	% frames as rows...
	% if n>m, disp([FUNCTION_NAME ': it is strongly suggested that m>n for [m,n] = size(M)']), end
	
	%%%%% set up initial guess for M
	Wx = state.W(1:2:end,:);
	Wy = state.W(2:2:end,:);
	if mod(m,2)==1 || ~isequal(Wx,Wy)
		if ~isequal(Wx,Wy), disp([FUNCTION_NAME ': suspicious - weight matrix does not have vertical pairing.']), end
		if mod(m,2)==1, disp([FUNCTION_NAME ': suspicious - M has an odd number of rows.']), end
		error([FUNCTION_NAME ':   Expecting vertically alternating x-y coords in M.'])
	end
	clear Wx Wy
	
	if ~isequal(state.B(:,r), ones(n,1))
        error([FUNCTION_NAME ': Expecting the last row of B to be entirely ones. ''SFM'' solutions only.'])
        %state.B(:,r) = 1;
    end
	
	% preallocate
	state.Hs = zeros(r-1,r-1,n);
	state.Gs = zeros(r-1,r-1,n,n);
	state.BIG_A = zeros((r-1)*(n-1),(r-1)*(n-1));
	state.BIG_B = zeros((r-1)*(n-1),1);
	state.xj_bar = zeros(m/2,r-1);
	state.t_bar = zeros(m,1);
	
	% image [row] weight sums
	state.sum_wjs = sum(state.W,2);
	state.r_sum_wjs = 1./state.sum_wjs;

    state.m = m;
    state.n = n;
    state.r = r;
    state.r_1 = r - 1;
    
    warning([FUNCTION_NAME ': the Brandt alternation scheme does not minimize the regularized error function; only sum-squared error is returned.']);
else
    % minimize
    
	% calculate points (for B)
	m_bar = sum(state.W.*state.M,2)./state.sum_wjs;
	
	for j=1:state.n
		tmp_sum = zeros(state.r_1,state.r_1);
		x = state.B(j,1:state.r_1)';
		for i=2:2:state.m
			P = state.A(i-1:i,1:state.r_1);
			tmp_sum = tmp_sum + (1-state.r_sum_wjs(i))*P'*P*state.W(i,j);
		end
		state.Hs(:,:,j) = tmp_sum;
		for jj=1:state.n
			tmp_sum = zeros(state.r_1,state.r_1);  
			for i=2:2:state.m
                tmp_sum =  tmp_sum + state.r_sum_wjs(i)*state.Hs(:,:,j)*state.W(i,jj); % this is the most expensive line of code
			end
			state.Gs(:,:,j,jj) = tmp_sum;
		end
	end
	
	for j=2:state.n
		j_indices = [1:state.r_1] + (j-2)*(state.r_1);
		for jj=2:state.n
			jj_indices = [1:state.r_1] + (jj-2)*(state.r_1);
			state.BIG_A(j_indices,jj_indices) = -state.Gs(:,:,j,jj) + state.Gs(:,:,j,1);
			if j==jj
                state.BIG_A(j_indices,jj_indices) = state.BIG_A(j_indices,jj_indices) + state.Hs(:,:,j);
			end
		end
	end
	
	for j=2:state.n
		j_indices = [1:state.r_1] + (j-2)*(state.r_1);
		tmp_sum = zeros(state.r_1,1);
		for i=2:2:state.m
			tmp_sum = tmp_sum + (1-state.r_sum_wjs(i))*state.A(i-1:i,1:state.r_1)'*(state.M(i-1:i,j)-m_bar(i-1:i))*state.W(i,j);
		end
		state.BIG_B(j_indices) = tmp_sum;
	end
	
	% check condition number of BIG_A
	if rcond(state.BIG_A)<1000*eps
		warning([FUNCTION_NAME ': BIG_A is nearly singular. Terminating now before inaccurate results take hold.'])
		return
	end
	
	state.BIG_X = state.BIG_A\state.BIG_B;
	
	tmp_sum = zeros(state.r_1,1);
	for j=2:state.n
		j_indices = [1:state.r_1] + (j-2)*(state.r_1);
		state.B(j,1:state.r_1) = state.BIG_X(j_indices,:)';
		tmp_sum = tmp_sum + state.B(j,1:state.r_1)';
	end
	state.B(1,1:state.r_1) = -tmp_sum';
	
	for i=2:2:state.m
        state.xj_bar(i/2,:) = sum((state.W(i,:)'*ones(1,state.r_1)).*state.B(:,1:state.r_1),1)/state.sum_wjs(i);
    end
	
	% calculate translations
	for i=2:2:state.m
        state.t_bar(i-1:i) = m_bar(i-1:i) - state.A(i-1:i,1:state.r_1)*state.xj_bar(i/2,:)'; 
    end
	state.A(:,state.r) = state.t_bar;
	
	% impute Mhat
	state.Mhat = amb_merge_matrices(state.W,state.M,state.A*state.B');
	
	% update cameras
	[U,S,V] = svd(state.Mhat);
	state.A(:,1:state.r_1) = U(:,1:state.r_1)*S(1:state.r_1,1:state.r_1);
	
	state.error = amb_fwmd_errfunc_sumsqrderr(state.A,state.B,state.M,state.W);
	    
end

