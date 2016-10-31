clear;
clc;

d = 100;
n = 100;                     %Data size
r = 4;                      %Rank

for kk = 1:100
% Generating data
    RU = randn(d,r);
    RV = randn(r,n);
    X_Ori = RU * RV;        %Original data
    Ind = randperm(d*n);
    X_Noi = X_Ori;
    
% Add 20% outliers into data
    p = floor(d*n*0.2);
    X_Noi(Ind(1:p)) = X_Noi(Ind(1:p)) + rand(1,p)* 10 - 5; %Add sparse noise
    
% Initialize U and V
    MedValX = median(abs(X_Noi(:)));
    MedValX = sqrt(MedValX/r);
    param.IniU = rand(d,r)*MedValX*2-MedValX;
    param.IniV = rand(n,r)*MedValX*2-MedValX;
    
% CWM implementation
    [U,V] = CWM(X_Noi,r,param);
    
% Record relative reconstruction errors (RRE)
    RRE(kk) = sum(sum(((X_Ori - U*V')).^2))/sum(sum(((X_Ori)).^2))
end

figure;
plot(RRE,'.');
xlabel('Experiments');
ylabel('RRE');
title('CWM performance');
