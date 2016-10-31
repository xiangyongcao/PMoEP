function [W TempU TempV] = maximizationW_PMoG(model,InW,InX,TempU,TempV,R,NumIter,param)
IND = find(InW(:)~=0);
k = size(R,2);
%mu = model.mu;
mu = zeros(1,k);
Sigma = model.Sigma;
w = model.Pi;
r = size(TempU,2);
W = zeros(size(InX));
C = zeros(size(InX));
for j = 1:k
    W(IND) = W(IND) + R(:,j)/(2*Sigma(j));
    C(IND) = C(IND) + R(:,j)*mu(j)/(2*Sigma(j));
end
C(IND) = C(IND)./W(IND);

if param.method == 1
    [TempU,TempV] = amb_mfwmd(TempU,TempV,InX-C,sqrt(W));
end

if param.method == 2
    [TempU,TempV] = EfficientMCL2(InX-C, sqrt(W), TempU,TempV, NumIter, 0.00000001);
end

if param.method == 3
    [M_est TempU TempV] = ALM(InX-C,sqrt(W),r,NumIter);
end

if param.method == 4
    [TempU,TempV,info,meanw]=weighted_pca(InX-C,sqrt(W),r-1,NumIter,0,TempU,TempV',mean(InX')');
    TempV = TempV';
    TempU(:,r) = meanw;
    TempV(:,r) = ones(size(TempV,1),1);
end