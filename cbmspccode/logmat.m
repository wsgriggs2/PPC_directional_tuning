function a = logmat(a)
[f v] = eig(a);
v = diag(log(diag(v)));
a = f*v*f';
