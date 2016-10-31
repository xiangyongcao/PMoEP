function y = logeppdf(x, lambda, p)
    y = log(p) + log(lambda)/p - log(2) - log(gamma(1/p))-lambda*(abs(x')).^(p);
