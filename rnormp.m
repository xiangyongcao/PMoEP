function y = rnormp(n,mu,sigmap,p,method)
% random number generator of the exponential power distribution
% input: n -  the number of the random number
%        mu - location parameter
%        sigmap - scale parameter
%        p - shape parameter
if nargin<5
    method = 'def';
end

if sigmap<=0
    error('sigmap must be positive');
end

if p>=1
    switch method
        case 'def'
            qg = gamrnd(1/p,p,1,n);
            z = qg.^(1/p);
            c1 = zeros(1,n); c2 = zeros(1,n);
            tmp = rand(1,n);
            c1(find(tmp<=0.5)) = -z(find(tmp<=0.5));
            c2(find(tmp>0.5)) = z(find(tmp>0.5));
            z = c1 + c2;
            y = -mu + z*sigmap;
            
        case 'chiodi'
            i = 0;
            y = zeros(1,n);
            while i<n
                u = 2*rand(1)-1; v = 2*rand(1)-1;
                z = abs(u)^p + abs(v)^(p/(p-1));
                if z<=1
                    i = i + 1;
                    y(i) = -u*(-p*log(z)/z)^(1/p);
                    y(i) = -mu + y(i)*sigmap;
                end
            end
    end
end

if p<1
   tau = (1/(p*sigmap^(p)))^(-1/p);
   w = 0.5*(1+p)*gamrnd(2+1/p,1,1,n) + 0.5*(1-p)*gamrnd(1+1/p,1,1,n);
   for i = 1:n
      c = tau*w(i)^(1/p);
      y(i) = sliceSampling(c, 100);
   end
end
% [nElems, centers] = hist(y, 500);
% bar(centers, nElems/n);axis([-10,60,0,0.01]);
