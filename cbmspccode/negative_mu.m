function [mu, grad_mu, hess_mu] = ...
    negative_mu(TransfMatrix,Sigma,Sigma_i,prob_i)

%NEGATIVE_MU  The negative of the mu-measure
%      
%[mu, grad_mu, hess_mu] = negative_mu(TransfMatrix,Sigma,Sigma_i,prob_i)
%
%TransfMatrix - (m x n) transformation matrix, where m and n are the
%dimensions of the feature space and original space, respectively
%
%Sigma - (n x n) overall covariances matrix
%
%Sigma_i - (n x n x C) class conditional covariance matrices, where C is
%the number of classes
%
%prob_i - (1 x C) prior probabilities of classes
%
%The mu-measure is described in:
%Z. Nenadic, Information discriminant analysis: Feature extraction with
%an information-theoretic objective, IEEE T. Pattern Anal., vol. 29 (8), 
%pp. 1394-1407, 2007.

%Author: Zoran Nenadic
%modified: 07/25/2007

%the number of classes
C = length(prob_i);

%the sizes of the feature and original spaces
[m n] = size(TransfMatrix);

%%GET THE FUNCTIONAL (NEGATIVE) MU
muo = 0;

for i = 1:C
  mu = muo + ...
       prob_i(i) * log(det(TransfMatrix * Sigma_i(:,:,i) * TransfMatrix'));
  muo = mu;
end

mu = 0.5 * [log(det(TransfMatrix * Sigma * TransfMatrix')) - muo];

%take the negative of mu
mu = (-1) * mu;

%%GET THE GRADIENT OF MU WITH RESPECT TO TRANSFORMATION MATRIX TransfMatrix
if nargout > 1
  
  grado = 0;
  
  for i = 1:C
    tmpmatr = TransfMatrix * Sigma_i(:,:,i) * TransfMatrix';
    grad_mu = grado + ...
              prob_i(i) * inv(tmpmatr) * TransfMatrix * Sigma_i(:,:,i);
    grado = grad_mu;
  end
  
  grad_mu = inv(TransfMatrix * Sigma * TransfMatrix') * ...
            TransfMatrix * Sigma - grado;
  
  %take the negative of the gradient
  grad_mu = (-1) * grad_mu;
end

%%GET THE HESSIAN OF MU WITH RESPECT TO TRANSFORMATION MATRIX TransfMatrix
if nargout > 2
  
  kk = 0;
  
  for i = 1:n
    for j = 1:m
      kk = kk + 1;
      b(kk) = kk;
      a(kk) = i + (j-1) * n;
    end
  end
  
  SparseMatrix = sparse(a,b,ones(1,m*n));
  
  hesso = 0;
  
  for i = 1:C
    temp1= TransfMatrix * Sigma_i(:,:,i);
    temp2 = inv(temp1 * TransfMatrix');
    Part_I = kron(eye(n),temp2) * kron(Sigma_i(:,:,i),eye(m))';
    Part_II = -kron(temp1',eye(m)) * kron(temp2',temp2) * ...
         [kron(temp1,eye(m)) + kron(eye(m),temp1) * SparseMatrix];
    hess_mu = hesso + prob_i(i) * [Part_I + Part_II];
    hesso = hess_mu;
  end
  
  temp1 = TransfMatrix * Sigma;
  temp2 = inv(temp1 * TransfMatrix');
  Part_I = kron(eye(n),temp2) * kron(Sigma,eye(m))';        
  Part_II = -kron(temp1',eye(m)) * kron(temp2',temp2) * ... 
       [kron(temp1,eye(m)) + kron(eye(m),temp1) * SparseMatrix];
  hess_mu = Part_I + Part_II - hesso;
  
  %take the negative of the hessian
  hess_mu = (-1) * hess_mu;
end
