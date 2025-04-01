function [Torthnorm, maxS] = orthonormalize(T)

%T - m x n matrix

[m n] = size(T);

%orthonormalize
[u s v] = svd(T);
is = 1./diag(s(1:m,1:m));
is = diag(is);
Torthnorm = is * u' * T;  
maxS = max(max(s));