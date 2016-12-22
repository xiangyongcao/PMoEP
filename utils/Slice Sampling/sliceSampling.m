function z = sliceSampling(c, maxIter)

z0 = 2 * c * (rand() - 0.5);
z = z0;

for iter = 1:maxIter
    u = (1 - abs(z/c)) * rand();
    z = 2*c*(1-u) * (rand() - 0.5);
end                                                                              

