function a = sqrtmat(a)
[f v] = eig(a);
v = diag(sqrt(diag(v)));
a = f*v*f';