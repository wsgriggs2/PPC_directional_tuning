function a = invsqrtmat(a)
[f v] = eig(a);
v = diag(1./sqrt(diag(v)));
a = f*v*f';
