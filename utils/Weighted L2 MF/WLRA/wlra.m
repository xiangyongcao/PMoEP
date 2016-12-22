function [X_t u s v] = wlra(W, X,IniX, K,Iter)
W = double(W);
X_t = IniX;
MaxIter = Iter;
k = 1;
while 1&k<MaxIter
    X_0 = X_t;
    [u s v] = lansvd(W.*X+(1-W).*X_t, K, 'L');
    X_t = u*s*v';
    if norm(X_t-X_0)/norm(X_0)<1e-4
        break;
    end
k = k+1;
end