function [M_est U_est V_est obj] = RobustApproximation_M_UV_TraceNormReg(M,W,r,maxIterIN,lambda,rho,signM)
%% Robust low-rank matrix approximation with missing data and outliers
% min |W.*(M-E)|_1 + lambda*|V|_*
% s.t., E = UV, U'*U = I
%
%Input:
%   M: m*n data matrix
%   W: m*n indicator matrix, with '1' means 'observed', and '0' 'missing'.
%   r: the rank of r
%   lambda: the weighting factor of the trace-norm regularization, 1e-3 in default.
%   rho: increasing ratio of the penalty parameter mu, usually 1.05.
%   maxInterIN: maximum iteration number of inner loop, usually 100.
%   signM: if M>= 0, then signM = 1, otherwise, signM = 0;
%Output:
%   M_est: m*n full matrix, such that M_est = U_est*V_est
%   U_est: m*r matrix
%   V_est: r*n matrix
%   L1_error: the L1-norm error of observed data only.
%% Normalization
scale = max(max(abs(M)));
%M = M/scale;

%% In-default parameters
[m n] = size(M); %matrix dimension
if nargin < 6
  rho = 1.05;
end
if nargin < 5
  lambda = 1e-3;
  signM = 0;
end
if nargin < 4
  maxIterIN = 100;
end
if nargin < 3
  disp('Please input the data matrix M, the indicator W and the rank r, and try again.');
end

maxIterOUT = 5000;
max_mu = 1e20;
mu = 1e-6;
M_norm = norm(M,'fro');
tol = 1e-8*M_norm;
tol = 1e-12;

cW = ones(size(W)) - W; %the complement of W.
display = 1; %display progress
%% Initializing optimization variables as zeros
E = randn(m,n);
U = randn(m,r);
V = randn(n,r);
Y = zeros(m,n); %lagrange multiplier
%% Start main outer loop
iter_OUT = 0;
objs=[];
while iter_OUT < maxIterOUT
  iter_OUT = iter_OUT + 1;
  if mod(iter_OUT,1000) == 0
      iter_OUT
  end
  
  itr_IN = 0;
  obj_pre = 1e20;
  %start inner loop
  while itr_IN < maxIterIN
    %update U
    
    U = (mu*E + Y)*V/(lambda*eye(r) + mu*V'*V);
    %  U = (E+Y/mu)*V/(V'*V+lambda/mu*eye(r));
    
    %update V
    V = (mu*E + Y)'*U/(lambda*eye(r) + mu*(U'*U));
    % V = (E+Y/mu)'*U/(U'*U+lambda/mu*eye(r));
    
    %update E
    UV = U*V';
    temp1 = UV - Y/mu;
    
    % l1
    %temp = M-temp1;
    %E = max(0,temp - 1/mu) + min(0,temp + 1/mu);
    
    %l1
    %E = (M-E).*W + temp1.*cW;
    E = 1/(mu+2)*(2*M+mu*temp1).*W + temp1.*cW;
    
    %evaluate current objective
    %obj_cur = sum(sum(abs(W.*(M-E)))) + lambda*sum(sigma0) + sum(sum(Y.*(E-UV))) + mu/2*norm(E-UV,'fro')^2;
    obj_cur = sum(sum(abs(W.*(M-E)).^2)) + lambda/2*(norm(U,'fro')^2 + norm(V,'fro')^2) + sum(sum(Y.*(E-UV))) + mu/2*norm(E-UV,'fro')^2;
    
    %check convergence of inner loop
    if abs(obj_cur - obj_pre) <= 1e-8*abs(obj_pre)
      break;
    else
      obj_pre = obj_cur;
      itr_IN = itr_IN + 1;
    end
  end
  
  leq = E - UV;
  stopC = norm(leq,'fro');
  if display
    % obj = sum(sum(abs(W.*(M-UV)))) + lambda*sum(sigma0);
    obj = sum(sum(abs(W.*(M-U*V')).^2)) + lambda/2*(norm(U,'fro')^2 + norm(V,'fro')^2);
    objs = [objs,obj];
  end
%   if display && (iter_OUT==1 || mod(iter_OUT,50)==0 || stopC<tol)
%     disp(['iter ' num2str(iter_OUT) ',mu=' num2str(mu,'%2.1e') ...
%       ',obj=' num2str(obj) ',stopALM=' num2str(stopC,'%2.3e')]);
%   end
  if stopC<tol
    break;
  else
    %update lagrage multiplier
    Y = Y + mu*leq;
    %update penalty parameter
    mu = min(max_mu,mu*rho);
  end
end

%% Denormalization
U_est = U;
V_est = V;

%U_est = sqrt(scale)*U; V_est = sqrt(scale)*V;
%L1_error = sum(sum(abs(W.*(scale*M-M_est))));
M_est = U_est*V_est';
obj = sum(sum(abs(W.*(M-M_est)).^2)) + lambda/2*(norm(U,'fro') + norm(V,'fro'));

end
