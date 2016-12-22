% M = amb_merge_matrices(W,M1,M2)
%
% Returns the 'alpha blend' of two matrices:
%
%      M  =  W.*M1  +  (1 - W).*M2

function M = amb_merge_matrices(W,M1,M2)

M = W.*M1 + (1-W).*M2;


