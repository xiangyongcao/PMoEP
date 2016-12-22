c = 1;
n = 100000;
maxIter = 100;

zz = zeros(1, n);

for i = 1:n
    zz(i) = sliceSampling(c, maxIter);
end

[nElems, centers] = hist(zz, 100);
bar(centers, nElems/n)

% [nElems, centers] = hist(E3(:), 100);
% bar(centers, nElems/numel(E3))