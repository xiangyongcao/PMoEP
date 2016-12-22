function [label,model,TW,OutU,OutV,llh,llh_BIC,p] = EM_PMoEP(InW,InX,r,param,p,lambda)
% Description: EM algorithm for fitting PMoEP model.
% USAGE: [label,model,TW,OutU,OutV,llh,llh_BIC,p] = EM_PMoEP(InW,InX,r,param,p,lambda)
%Input:
%   InW: d x n x param.k indicator matrices
%   InX: d x n input data matrix
%   r:   the rank
%   param:
%      --param.maxiter: maximal iteration number
%      --param.OriX: ground truth matrix
%      --param.InU: initialized factorized matrice U
%      --param.InV: initialized factorized matrice V
%      --param.k: the number of mixture components
%      --param.display: display the iterative process
%      --param.tol: the thresholding for stop
%   p: the candidates components
%   lambda: the tuning parameter

%Output:
%   label: the labels of the noises
%   model: model.eta, the precisions of the different EPs
%          model.Pi,the mixing coefficients
%   W: d x n weighted matrix
%   OutU: the final factorized matrix U
%   OutV: the final factorized matrix V
%   llh:  the log likelihood
%   llh_BIC:  the log likelihood used in BIC criterion
%   p: the selected components

% Author: Xiangyong Cao(caoxiangyong45@gmail.com)
% If you have any quesion, please contact Xiangyong Cao

% Reference paper: 
% @InProceedings{Cao_2015_ICCV,
% author = {Cao, Xiangyong and Chen, Yang and Zhao, Qian and Meng, Deyu and Wang, Yao and Wang, Dong and Xu, Zongben},
% title = {Low-Rank Matrix Factorization Under General Mixture Noise Distributions},
% journal = {The IEEE International Conference on Computer Vision (ICCV)},
% month = {June},
% year = {2015}
% }

% version 1.0, date: 12-24-2015

% initialization
[d,n] = size(InX);
if (~isfield(param,'maxiter'))
    maxiter = 100;
else
    maxiter = param.maxiter;
end

if (~isfield(param,'OriX'))
    OriX = InX;
else
    OriX = param.OriX;
    clear param.OriX;
end

IND = find(InW(:) ~= 0);
if (~isfield(param,'InU'))
    s = median(abs(InX(IND)));
    s = sqrt(s/r);
    if min(InX(IND)) >= 0
        InU = rand(d,r)*s;
    else
        InU = rand(d,r)*s*2-s;
    end
else
    InU = param.InU;
end

if (~isfield(param,'InV'))
    if min(InX(IND)) >= 0
        InV = rand(n,r)*s;
    else
        InV = rand(n,r)*s*2-s;
    end
else
    InV = param.InV;
end

if (~isfield(param,'k'))
    k = 3;
else
    k = param.k;
end

if (~isfield(param,'display'))
    display = 0;
else
    display = param.display;
end

if (~isfield(param,'tol'))
    tol = 1e-7;
else
    tol = param.tol;
end

param.method = 2;
%Initialize the parameters of MoEP
IND = find(InW(:) ~= 0);
tempX=InX(IND);

if length(p)~=length(find(p==2))
    R = initialization_PMoEP(tempX',k,'random'); % initialize label
    [~,label(1,:)] = max(R,[],2);
    R = R(:,unique(label));
    eta = 10*rand(1,k); % eta = 10*ones(1,k);
    nk = sum(R,1);
    Pi =  nk/size(R,1);
    model.eta = eta;
    model.Pi = Pi;
else
    R = initialization_PMoG(tempX',k);
    [~,label(1,:)] = max(R,[],2);
    R = R(:,unique(label));
    model.Sigma = rand(1,k);
    nk = sum(R,1);
    model.Pi = nk/size(R,1);
    model.mu = zeros(1,k);
end

% llh = -inf(1,maxiter);
converged = false;
TempU = InU;
TempV = InV;
TempX = TempU * TempV';
Error = InX(:)-TempX(:);
Error = Error(IND);

t = 1;

%Initialized E Step
if length(p)~=length(find(p==2))
    [R, llh(t),llh_BIC] = expectation(Error',model, p, lambda);
else
    [R, llh(t),llh_BIC] = expectation_PMoG(Error',model,lambda);
end

while ~converged && t < maxiter
   
    t = t+1;
   
    % M Step
    if length(p)~=length(find(p==2))
        [model,~,p,echo] = maximizationModel(Error',R,model.eta,p,lambda);
    else
        [model,R,echo,k] = maximizationModel_PMoG(Error',R,lambda);
        p = 2*ones(1,k);
    end
    
    if echo==1
        OutU = InU; OutV = InV; llh=llh(t-1);
        llh_BIC=llh; label=[]; W = InW;
        return;
    end
    
    % E Step
    if length(p)~=length(find(p==2))
        [R, llh(t),~] = expectation(Error',model, p, lambda);
        L1 = llh(t);
    else
        Sigma = 1./(2*model.eta);
        model.Sigma = Sigma;
        [R, llh(t),llh_BIC] = expectation_PMoG(Error',model,lambda);
        L1 = llh(t);
    end
    
    % M Step
    if length(p)~=length(find(p==2))
        [TempU, TempV] = ALM_UV(InW,model,InX,r,R,p);
    else
        [TW TempU TempV] = maximizationW_PMoG(model,InW,InX,TempU,TempV,R,100,param);
    end
    TempX = TempU * TempV';
    Error = InX(:) - TempX(:);
    Error = Error(IND);
    
    % E Step
    if length(p)~=length(find(p==2))
        [R, llh(t), llh_BIC] = expectation(Error',model, p, lambda);
        L2 = llh(t);
    else
        [R, llh(t),llh_BIC] = expectation_PMoG(Error',model,lambda);
        L2 = llh(t);
    end
    
    [~,label] = max(R,[],2);
    u = unique(label);   % non-empty components
    if size(R,2) ~= size(u,2)
        R = R(:,u);   % remove empty components
        R = R./repmat(sum(R,2),1,size(R,2));
        if length(p)~=length(find(p==2))
            model.eta = model.eta(u);
            model.Pi = model.Pi(u);
            model.Pi = model.Pi/sum(model.Pi);
            p = p(u);
        else
            model.mu = model.mu(u);
            model.Sigma = model.Sigma(u);
            model.Pi = model.Pi(u);
        end
    end
    converged = abs(llh(t)-llh(t-1)) < tol;
    k = length(u);
    if mod(t,10)==0
    disp(['Iteration ',num2str(t),': there are ',num2str(k),...
        ' EP components and the corresponding p are ',num2str(p),'. Likelihood in this step is ',num2str(L2),'.']);
    end
end

if length(p)~=length(find(p==2))
    TW = zeros(d,n,k);
    for l = 1:k
        Temp = zeros(d,n);
        Temp(IND)=R(:,l);
        TW(:,:,l) = (model.eta(l)*Temp).^(1/p(l));
    end
end

[~,label] = max(R,[],2);
OutU = TempU; OutV = TempV;

if ~display
    disp(['The likelihood in this step is ',num2str(L1),' and ',num2str(L2),';']);
    if length(p)~=length(find(p==2))
        disp(['There are ',num2str(k),' hylaplace noises mixed in data']);
        disp(['The selected p is ',num2str(p)]);
        disp(['with precisions ',num2str(model.eta)]);
        disp(['with weights ',num2str(model.Pi),'.']);
    else
        disp(['There are ',num2str(k),' Gaussian noises mixed in data']);
        disp(['with means ',num2str(model.mu)]);
        disp(['with variances ',num2str(model.Sigma)]);
        disp(['with weights ',num2str(model.Pi),'.']);
    end
    if (~isfield(param,'orimissing'))
        disp(['Relative reconstruction error ', num2str(sum(sum(((OriX - TempU*TempV')).^2))/sum(sum((OriX).^2)))]);
    else
        disp(['Relative reconstruction error ', num2str(sum(sum((InW.*(OriX - TempU*TempV')).^2))/sum(sum((InW.*OriX).^2)))]);
    end
    disp(['L2 RMSE is ', num2str(sqrt(mean(Error.^2)))]);
    disp(['L1 RMSE is ', num2str(mean(abs(Error)))]);
end
if converged
    fprintf('Converged in %d steps.\n',t-1);
else
    fprintf('Not converged in %d steps.\n',maxiter);
end
