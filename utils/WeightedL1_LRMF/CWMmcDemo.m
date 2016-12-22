clear;
clc;

d = 100;
n = 100;                    %Data size
r = 4;                      %Rank

for kk = 1:100
% Generating data
    RU = randn(d,r);
    RV = randn(r,n);
    X_Ori = RU * RV;        %Original data
    Ind = randperm(d*n);
    p1 = floor(d*n*0.2);
    W = ones(d,n);
    W(Ind(1:p1)) = 0;       %Indicator matrix
    X_Noi = X_Ori;
    
% Add 20% missing components into data
    X_Noi = W.*X_Noi; 
    
% Add 20% outliers into data
    p2 = floor(d*n*0.2);
    X_Noi(Ind(p1+1:p1+p2)) = X_Noi(Ind(p1+1:p1+p2)) + rand(1,p2)* 10 - 5;
    
% Initialize U and V
    MedValX = median(abs(X_Noi(Ind(p1+1:end))));
    MedValX = sqrt(MedValX/r);
    param.IniU = rand(d,r)*MedValX*2-MedValX;
    param.IniV = rand(n,r)*MedValX*2-MedValX;
    
% CWM implementation
    [U,V] = CWMmc(X_Noi,W,r,param);
    
% Record relative reconstruction errors (RRE)
    RRE(kk) = sum(sum(((X_Ori - U*V')).^2))/sum(sum(((X_Ori)).^2))
end

figure;
plot(RRE,'.');
xlabel('Experiments');
ylabel('RRE');
title('CWM performance');
