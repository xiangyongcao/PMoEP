function [model,R,echo,k] = maximizationModel_PMoG(X,R,lambda)

[d,n] = size(X);
k = size(R,2);

% nk = sum(R,1);
mu = zeros(1,k);%fix mu to zero
Sigma = rand(1,k);

% mu = bsxfun(@times, X*R, 1./nk);

Tmp1 = 1/(1-k*lambda*2)*(sum(R)/n-lambda*2);
w = max(zeros(1,k),Tmp1);
if length(find(w==0))~=length(w)
    w = w/sum(w);
end

ind = find(w==0); len_ind = length(ind);

if k > len_ind
    
    % delete component
    k = k - len_ind;
    w(ind) = []; R(:,ind) = []; mu(:,ind)=[]; Sigma(:,ind) = [];
    
    % Dealing with the case a row of R is all 0 after the delete
    index = find(sum(R,2)==0); len = length(index);
    
    % clustering the noise again corresponding the posterior is all 0
    Error_index = X(index);
    model.mu = mu; model.Sigma = Sigma; model.Pi = w;
    [R_index] = expectation_PMoG(Error_index',model,lambda);
    R(index,:) = R_index;
    R = R./repmat(sum(R,2),1,k);        % if sum(isnan(R))~=0 pause;end
    
    sqrtR = sqrt(R); nk = sum(R,1);
     
    for i = 1:k
        Xo = bsxfun(@minus,X,mu(i));
        Xo = bsxfun(@times,Xo,sqrtR(:,i)');
        Sigma(i) = (Xo*Xo')/nk(i);
        Sigma(i) = Sigma(i)+(1e-6); % add a prior for numerical stability
    end
    model.mu = mu; model.Sigma = Sigma; model.eta = 1./(2*Sigma); model.Pi = w; echo = 0;
    
else
    model.mu = mu; model.Sigma = Sigma; model.eta = 1./(2*Sigma); model.Pi = w; echo = 1;
end