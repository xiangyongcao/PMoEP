Usage of damped_wiberg:

This MATLAB function factorizes a given matrix Y into U*V', a product
of two smaller matrices of rank r. To be specific, it minimizes ||H
\odot (Y - U*V')||_2 for specified data Y, H, and r, where Y is the
matrix to be factorized, H is the binary matrix indicating existing
(1) and missing (0) components of Y, and r is the number of columns of
U and V. For more detail, see the source code. 

For the details of the implementation, please see the following
reference. Please cite it in your work when you use our code. 

[1] Takayuki Okatani, Takahiro Yoshida, Koichiro Deguchi: Efficient
algorithm for low-rank matrix factorization with missing components
and performance comparison of latest algorithms. Proc, ICCV 2011:
842-849.

Example)

1. Load your data. (These data are from
http://www.robots.ox.ac.uk/~abm/. Please download the original when
you use them for your study.)

>> load('M.txt');
>> load('mask.txt');

2. Factorize the matrix.

>> [U,V]=damped_wiberg(M, mask, 4);

Done. The algorithm should ALWAYS converge to the global minimum for
about 50-300 iterations. The results are such that M <-> U*V'.
