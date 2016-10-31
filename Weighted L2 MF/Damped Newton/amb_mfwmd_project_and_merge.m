function state = amb_mfwmd_project_and_merge(state,r)

% state = amb_mfwmd_project_and_merge(state[,r])
%
% Perform one iteration of the project and merge algorithm. Provide all
% arguments for initialization of 'state' struct. The fields that are
% required to already exist are: 
%
%   A, B, M, W
%
% Citation: Guerreiro and Aguiar, "3D Structure from Video Streams with
% Partially Overlapping Images", ICIP'02
%
% See 'amb_mfwmd'

FUNCTION_NAME = 'amb_mfwmd_project_and_merge';

if nargin==2
    % initialize

    state.fstr = 1:r; % first 'r' columns
    % intial Mhat
    % (imputing the missing elements with an average value can cause initial increase in error)
    state.Mhat = amb_merge_matrices(state.W,state.M,state.A*state.B');
    
    warning([FUNCTION_NAME ': the project and merge scheme does not minimize the regularized error function; only sum-squared error is returned.']);
else
    % minimize
    
	% Project
	[U,S,V] = svd(state.Mhat);
	state.A = U(:,state.fstr)*S(state.fstr,state.fstr);
	state.B = V(:,state.fstr);
	
	% Merge
	state.Mhat = amb_merge_matrices(state.W,state.M,state.A*state.B');

    state.error = amb_fwmd_errfunc_sumsqrderr(state.A,state.B,state.M,state.W);
end
