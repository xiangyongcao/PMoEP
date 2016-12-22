function [R, llh, llh_BIC] = expectation(x, model, p, C)
eta = model.eta;  
Pi = model.Pi;

n = length(x);
k = length(Pi);
logRho = zeros(n,k);
epsilon = 1e-6;

for l = 1:k
    logRho(:,l) = logeppdf(x,eta(l),p(l));
end
logRho = bsxfun(@plus,logRho,log(Pi));
T = logsumexp(logRho,2);
llh_BIC = sum(T);
llh = sum(T)-n*C*2*(sum(log(epsilon+Pi))-k*log(epsilon)); % uncomplete_penalized_loglikelihood
% llh = sum(T)/n - n*C*(sum(log(epsilon+Pi))-k*log(epsilon)); % loglikelihood %%divide n is not correct
logR = bsxfun(@minus,logRho,T);
R = exp(logR);