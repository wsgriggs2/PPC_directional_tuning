function [mu, grad_mu] = ...
    negative_mu_cg(TransfMatrix,Parameters)

%NEGATIVE_MU_CG  The negative of the mu-measure modified for
%compatibility with conj_grad.m (written by Hans Bruun Nielsen) and
%downloadable from http://www2.imm.dtu.dk/~hbn/Software/
%      
%[mu, grad_mu] = negative_mu_cg(TransfMatrix,Sigma,Sigma_i,prob_i)
%
%TransfMatrix - (1 x (m * n)) transformation matrix, where m and n are the
%dimensions of the feature space and original space, respectively.
%Note TransfMatrix has to be reshaped into a vector form to accommodate
%the constraints of conj_grad.m
%
%Parameters - (1 x 1) structure:
%
%Parameters.Sigma - (n x n) overall covariances matrix
%
%Parameters.Sigma_i - (n x n x C) class conditional covariance matrices, 
%where C is the number of classes
%
%Parameters.prob_i - (1 x C) prior probabilities of classes
%
%The mu-measure is described in:
%Z. Nenadic, Information discriminant analysis: Feature extraction with
%an information-theoretic objective, IEEE T. Pattern Anal., vol. 29 (8), 
%pp. 1394-1407, 2007.

%Author: Zoran Nenadic
%modified: 09/26/2007

Sigma = Parameters.Sigma;
Sigma_i = Parameters.Sigma_i;
prob_i = Parameters.prob_i;

%the size of data vector
n = size(Sigma,1);

%the size of feature space
m = length(TransfMatrix)/n;

%reshape back the matrix to take advantage of matrix manipulations
TransfMatrix = reshape(TransfMatrix,m,n);

%the number of classes
C = length(prob_i);

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
grad_mu = (-1) * grad_mu;   %(m x n) matrix

%reshape gradient to the vector form
grad_mu = reshape(grad_mu,1,m*n);


