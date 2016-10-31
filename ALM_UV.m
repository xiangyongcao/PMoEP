function [U,V] = ALM_UV(W,model,InX,r,R,p)

[d,n] = size(InX);

% Initialize U and V
MedValX = median(abs(InX(:)));
MedValX = sqrt(MedValX/r);
para.IniU = rand(d,r)*MedValX*2-MedValX;
para.IniV = rand(n,r)*MedValX*2-MedValX;
para.display = 0;
para.maxiter = 10;
[U,V] = CWM(InX,r,para);

% ALM initialization
E = InX - U*V';
Y = zeros(d,n);
rho = 2; rho0 = 1;
alpha = 1.2; rho_max = 5*rho0; %1.25
 IND=find(W~=0);
 
while rho<rho_max
    % update U,V
    temp= InX-E- Y/rho;
    [U,V,~] = EfficientMCL2(temp, W, U,V,10,1e-6);
    
    % update E
    Tau = model.eta./rho;
    y = InX-U*V'- Y/rho; 
    E(IND) = Newton_lp(y(IND),Tau,R,p,0.001);
    
    % update Y
    G=E + U*V'-InX;
    Y(IND) = Y(IND) + rho*G(IND);
    
    % update rho
    rho = alpha*rho;
end

function x=Newton_lp(y,lamda,R,p,eps)
% solve: min_x 0.5*£¨y-x£©^2+sigma{lamda(i).*R(:,i).*abs(x)^p(i)}
% If y>0, result is x. 
% If y<0, result is -x.
t=y;
y=abs(y);
x0=y./(max(lamda)+1);
x1=x0-fun(x0,y,lamda,R,p)./gfun(x0,y,lamda,R,p); 
k=0;
while k<20 & norm(fun(x1,y,lamda,R,p))>eps 
   alpha=1;
   x0=x1;
   x1=x0-alpha*fun(x0,y,lamda,R,p)./gfun(x0,y,lamda,R,p);
    while abs(fun(x0,y,lamda,R,p))<abs(fun(x1,y,lamda,R,p)) %adjust step size
    alpha=alpha*0.8;
    x1=x0-alpha*fun(x0,y,lamda,R,p)./gfun(x0,y,lamda,R,p);
    end
    k=k+1;
end

x=zeros(size(x1));
L1=find(x1>0);
L2=find(yfun(x1,y,lamda,R,p)<yfun(x,y,lamda,R,p));
L=intersect(L1,L2);
x(L)=x1(L);
x=x.*sign(t);

% define objective function and first order derivative
function f=fun(x,y,lamda,R,p)
k=length(lamda);
f=x-y;
for i=1:k
    f=f+lamda(i).*R(:,i).*p(i).*abs(x).^(p(i)-1).*sign(x);
end

function g=gfun(x,y,lamda,R,p)
k=length(lamda);
g=ones(size(y));
for i=1:k
    g=g+lamda(i).*R(:,i).*p(i).*(p(i)-1).*abs(x).^(p(i)-2);
end

% original function
function q=yfun(x,y,lamda,R,p)
k=length(lamda);
q=0.5*(x-y).^2;
for i=1:k
    q=q+lamda(i).*R(:,i).*abs(x).^p(i);
end
