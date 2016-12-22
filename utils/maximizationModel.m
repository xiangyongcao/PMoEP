function [model,R,p,echo] = maximizationModel(x,R,eta,p,C)

n = length(x);
k = size(R,2);

% update Pi
Tmp1 = 1/(1-2*k*C)*(sum(R)/(n)-2*C);
Pi = max(zeros(k,1),Tmp1')';

if length(find(Pi==0))~=length(Pi)
    Pi = Pi/sum(Pi);
end

ind = find(Pi==0); len_ind = length(ind);

if k > len_ind
    % delete component
    ind = find(Pi==0);
    k = k - length(ind);
    Pi(ind) = [];
    R(:,ind) = [];
    R = R./repmat(sum(R,2),1,k);
    eta(ind) = [];
    p(ind) = [];
    
    N = sum(R);
    
    % update lambda
    for l = 1:k
        tmp = sum(R(:,l)'.*((abs(x)).^(p(l))));
        eta(l) = N(l)/(p(l)*tmp);
    end
    
    model.eta = eta;
    model.Pi = Pi;
    echo = 0;
else
    model.eta = eta;
    model.Pi = Pi;
    echo = 1;
end