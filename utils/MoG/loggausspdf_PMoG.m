function y = loggausspdf_PMoG(X, mu, Sigma)
d = size(X,1);
X = bsxfun(@minus,X,mu);
    
% if (Sigma==[])|(isnan(Sigma)==1)
%     Sigma = 1e-6;
% end

[U,p]= chol(Sigma);

% if p ~= 0
%     error('ERROR: Sigma is not PD.');
% end

% disp(['X is', num2str(X)]);
% disp(['U is', num2str(U)]);
% disp(['Sigma is', num2str(Sigma)]);
% disp(['p is', num2str(p)]);
 
Q = U'\X;
q = dot(Q,Q,1);  % quadratic term (M distance)
c = d*log(2*pi)+2*sum(log(diag(U)));   % normalization constant
y = -(c+q)/2;

% y1 = -0.5*(log(2*pi)+log(Sigma));
% y2 = -0.5*(X-mu).^2/Sigma;
% y = y1+y2;