function [R, llh, llh_BIC] = expectation_PMoG(X, model,lambda)
% mu = model.mu;
Sigma = model.Sigma;
w = model.Pi;


n = length(X);
k = length(w);
logRho = zeros(n,k);
epsilon = 1e-6;

for l = 1:k
    logRho(:,l) = loggausspdf_PMoG(X,0,Sigma(l))';
end
logRho = bsxfun(@plus,logRho,log(w));
T = logsumexp(logRho,2);
llh_BIC = sum(T);
llh = sum(T)-n*lambda*2*(sum(log(epsilon+w))-k*log(epsilon)); % uncomplete_penalized_loglikelihood
logR = bsxfun(@minus,logRho,T);
R = exp(logR);