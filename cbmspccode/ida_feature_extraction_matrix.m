function [Tida, Muida] = ...
    ida_feature_extraction_matrix(m,Train,Group,Method, ...
                                  Tol,MaxIter,InitCondition,Nruns)
  
%IDA_FEATURE_EXTRACTION_MATRIX  Optimization scheme for finding IDA
%feature extraction matrix. The details of the IDA method can be found
%in:
%Z. Nenadic, Information discriminant analysis: Feature extraction with
%an information-theoretic objective, IEEE T. Pattern Anal., vol. 29 (8), 
%pp. 1394-1407, 2007.
%      
%[Tida, Muida] = ida_feature_extraction_matrix(m,Train,Group,Method, ...
%                                Tol,MaxIter,InitCondition,Nruns)
%
%Input arguments:
%---------------
%
%m - (1 x 1) the size of feature space (chosen by the user)
%
%Train - (Nt x n) is the statistical data of interest, where Nt is the
%number of samples (or training instances) in the data set, and n is the
%size of the data vector (number of variables). To avoid overfitting, 
%feature extraction must be based on the training data, hence the
%name "Train".
%
%Group - (Nt x 1) is a vector of class labels corresponding to the rows
%in Train. The convention is as follows: first class is labeled as 0,
%second class is labeled as 1, ..., C-th class is labeled as (C-1), 
%where C is the number of classes   
%
%Method - is a string: 'tr' for the trust-region method
%                      'cg' for the conjugate gradient method
%
%Tol - (1 x 2) vector of tolerances: Tol(1) tolerance for the fcn. value
%                                    Tol(2) tolerance for the fcn. argument
%
%MaxIter - (1 x 1) maximum number of iterations
%
%InitCondition - is a string: 'random' random initial feature extraction matrix
%                             'lda' linear discriminant analysis matrix
%                             'che' the matrix proposed by Loog & Duin in
%                             IEEE TPAMI, 26, 732-739, 2004. 
%
%Nruns - (1 x 1) the number of optimization runs. Choosing large Nruns
%will slow down the process, whereas leaving Nruns small may result in a
%local optimum. Local solution from previous run is randomly perturbed
%and used as an initial condition for the next run.
%
%Output arguments:
%----------------
%
%Tida - (m x n) is the feature extraction matrix
%
%Muida - (1 x 1) is the value of the class-separability function called
%the "mu measure"
%  
%[Tida, Muida] = ida_feature_extraction_matrix(m,Train,Group,Method, ...
%                                  Tol,MaxIter)
%or
%
%[Tida, Muida] = ida_feature_extraction_matrix(m,Train,Group,Method, ...
%                                  Tol,MaxIter,[])
%or
%
%[Tida, Muida] = ida_feature_extraction_matrix(m,Train,Group,Method, ...
%                                  Tol,MaxIter[],[])
%
%will call the function with the default parameters T0 (random matrix)
%and Nruns (1 run) 
  
%Author: Zoran Nenadic
%modified: 09/18/2007

%check the number of arguments
if nargin < 6
  error('arguments missing; type help ida_feature_extraction_matrix');
elseif nargin == 6  
  InitCondition = 'random';
  Nruns = 1;
elseif nargin == 7
   Nruns = 1;
elseif nargin > 8
  error('too many arguments; type help ida_feature_extraction_matrix');
end 

%the size of the data space  
[Nt n] = size(Train);

%the number of classes
C = length(unique(Group));

%calculate the sample statistics 
Mean = zeros(n,1);
Sigma_w = zeros(n,n);
Sigma_b = zeros(n,n);

for i = 1:C
    
    ind = find(Group == i-1);
    
    prob_i(i) = length(ind)/Nt;
    Mean_i(:,i) = mean(Train(ind,:),1)';
    Sigma_i(:,:,i) = cov(Train(ind,:));
    
    %overall mean
    Mean = Mean + prob_i(i) * Mean_i(:,i);
    
    %within class matrix
    Sigma_w = Sigma_w + prob_i(i) * Sigma_i(:,:,i);

    %between class matrix
    Sigma_b = Sigma_b + prob_i(i) * [Mean_i(:,i) * Mean_i(:,i)'];

end

Sigma_b = Sigma_b - Mean * Mean';
Sigma = Sigma_w + Sigma_b;

switch num2str(InitCondition)
 case 'random'
  T0 = randn(m,n);
 case 'lda'
  T0 = linear_disc_analysis(Sigma_w,Sigma_b,m);
 case 'che'
  T0 = ...
      acc_method(Sigma_w, Sigma_b, Sigma_i, Mean_i, prob_i, m, C);
 otherwise
  %default is a randomized T0
  T0 = randn(m,n);
end

if isempty(Nruns)
  Nruns = 1;
end

%trust-region method
if strcmp(Method,'tr')
  disp('*** running trust-region method ***')
  OPT = optimset('fminunc');
  OPT = optimset(OPT,'GradObj','on','Hessian','on','Display','off', ...
                 'TolFun',Tol(1),'TolX',Tol(2),'MaxIter',2000);
  
  Tbest = T0;   %log the best result
  muo = (-1) * negative_mu(T0,Sigma,Sigma_i,prob_i);
  Tinit = T0;   %set initial condition
  
  for j = 1:Nruns
    
    fprintf('Run: %d\n',j)
    
    %orthonormalize T0
    [Tinit, maxS] = orthonormalize(Tinit);
    
    [Topt, mu_opt] = fminunc(@(T) negative_mu(T,Sigma,Sigma_i,prob_i), ...
                             Tinit,OPT);

    if abs(mu_opt) > muo
        fprintf('better solution found: %1.13f\n',(-1)*mu_opt);
        muo = (-1) * mu_opt;
        Tbest = Topt;   %update the best result
    end

    %random restart
    Tinit = Tbest + randn(m,n) * maxS;
    
  end
  fprintf('final solution: %1.13f\n',muo);
  %fprintf('\n')

  %orthonormalize best result
  Tida = orthonormalize(Tbest);
  Muida = muo;

%conjugate gradient method  
elseif strcmp(Method,'cg')
  disp('*** running conjugate gradient method ***')
  %create a structure (constraint imposed by conj_grad.m)
  Parameters.Sigma = Sigma;
  Parameters.Sigma_i = Sigma_i;
  Parameters.prob_i = prob_i;
  
  %set options for conjugate gradient method (see conj_grad.m)
  OPT(1) = 1; %1 - Fletcher-Reeves, 2 - Polak-Ribiere
  OPT(2) = 1; %1 - exact line search, 2 - soft line search
  OPT(3) = 1; %upper bound on initial step
  [muo go] = negative_mu_cg(reshape(T0,1,m*n),Parameters); 
  OPT(4) = Tol(1) * norm(go);  
  OPT(5) = Tol(2); 
  OPT(6) = MaxIter;
  
  %set options for line search (see linesearch.m)
  OPT(7) = 0.99;   %stopping criterion 
  OPT(8) = 1e-3;   %stopping criterion
  OPT(9) = 5;     %maximum number of iterations
  
  Tbest = T0;   %log the best result
  muo = (-1) * muo;
  Tinit = T0;   %set initial condition
  
  for j = 1:Nruns
    fprintf('Run: %d\n',j)
    
    %orthonormalize T0
    [Tinit, maxS] = orthonormalize(Tinit);
    
    %reshapte the matrix to make a vector
    Tinit = reshape(Tinit,1,m*n);
    
    [Topt, Info] = ...
        conj_grad('negative_mu_cg',Parameters,Tinit,OPT);
    
    if abs(Info(1)) > muo
      muo = (-1) * Info(1);
      fprintf('better solution found: %1.13f\n',muo);
      
      Tbest = Topt;   %update the best result
    end

    %random restart
    Tinit = Tbest + randn(1,m*n) * maxS;
    
  end
  fprintf('final solution: %1.13f\n',muo);
  
  %orthonormalize best result
  Tida = orthonormalize(reshape(Tbest,m,n));
  Muida = muo;
else
  error('unknown optimization method')
end


%-------- SUBROUTINES ----------------

%-------- LDA ------------------------
function [Tlda] = linear_disc_analysis(Sigma, Sigma_b, m)
  
Rsb = rank(Sigma_b);
Rs = rank(Sigma);
n = size(Sigma,1);

if m > Rsb | n > Rs 
  disp('WARNING: LDA not feasible, switching to random matrix');
  Tlda = randn(m,n);
else
  [V D] = eig(inv(Sigma) * Sigma_b);
  D = diag(D);
  [D Ind] = sort(D);
  V = V(:,Ind);
  V = fliplr(V);
  Tlda = V(:,1:m)';
end
 
%-------- CHE ------------------------
function [Tche] = ...
      acc_method(Sigma_w, Sigma_b, Sigma_i, Mean_i, prob_i, m, C)
   
n = size(Sigma_w,1);
  
%sphering
rootSw = sqrtmat(Sigma_w);   %transforms the variables so Sigma_w = I;
W = inv(rootSw);             %whitening matrix so that Sigma_w = I;

for i = 1:C
    Zi(:,:,i) = W * Sigma_i(:,:,i) * W;  %transformed si's
end 

Sche = zeros(n,n);

for i = 1:C-1
    for j = i+1:C
        Mij = Mean_i(:,i) - Mean_i(:,j);
        pi_i = prob_i(i)/(prob_i(i) + prob_i(j));
        pi_j = prob_i(j)/(prob_i(i) + prob_i(j));
        Zij = W * [pi_i * Sigma_i(:,:,i) + pi_j * Sigma_i(:,:,j)] * W; 
        Sche = Sche + prob_i(i) * prob_i(j) * rootSw * [ ...
            (inv(sqrtmat(Zij))) * W * Mij * ...
            Mij' * W * (inv(sqrtmat(Zij))) + ...
            1/(pi_i*pi_j)*[logmat(Zij) - pi_i * logmat(Zi(:,:,i)) - ...
            pi_j * logmat(Zi(:,:,j))]] * rootSw;
    end
end

[V D]= eig(inv(Sigma_w)*Sche);
[D iD] = sort(diag(D));
V = fliplr(V(:,iD));
V = V(:,1:m);
Tche = V/norm(V);
Tche = Tche';

function R = sqrtmat(M)
%square root of a matrix
%note that root is a symmetric matrix

[V D]= eig(M);
R = V * D.^0.5 * V'; 

function R = logmat(M)
%logarithm of a matrix

[V D]= eig(M);
d = diag(D);
logd = log(d);
D = diag(logd);
R = V * D * V';