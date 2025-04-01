function [OptSubInd, OptClass] = choose_subspace...
    (TestData,TrainData,TrainLabels,Subspace,Prior)

%Input arguments:
%1) TestData (1 x N) vector (a single sample of test data), with
%N being the dimension of the data (number of attributes)
%2) TrainData (n x N) matrix, where n is the number of samples (instances),
%and N is the dimension of the data.
%3) TrainLabels (n x 1) vector of class labels corresponding to TrainData
%4) Subspace (1 x C) cell, whrere C is the number of classes. These are the
%subspaces returned by dataproc_func_cpca.m
%5) Prior (string) specifies the way prior (class) probabilities are
%calculated in the method. Choices: 'empirical' (priors match the,relative
%frequencies of classes), or 'uniform' (priors are uniform-default)
%
%Output arguments:
%1) OptSubInd (1 x 1) the index of optimal subspace
%2) OptClass (1 x 1) is the label of the optimal class (this is not 
%recommended as a classification tool).  

%Zoran Nenadic DSc
%06/21/2010

classes = unique(TrainLabels);
Nclass = length(classes);
Nobs = length(TrainLabels);
for c = 1:Nclass
    NtrialA(c) = length(find(TrainLabels == classes(c)));
end

switch Prior
    case {'empirical'}
        P_i = NtrialA./Nobs;    %prior class probabilities
    otherwise
        P_i = ones(1,Nclass);     %default to uniform prior
        P_i = P_i/sum(P_i);
end
%disp('---------------')
for i = 1:Nclass
    %disp(['subspace #' num2str(i)])
    TestFeature = TestData * Subspace{i};
    TrainFeature = TrainData * Subspace{i};
    Posterior = zeros(1,Nclass);
    for j = 1:Nclass
        indc = find(TrainLabels==classes(j));
        mu = mean(TrainFeature(indc,:),1);
        Sigma = cov(TrainFeature(indc,:));
        TF = TestFeature - mu;
        Posterior(j)=-0.5*TF*inv(Sigma)*TF' - 0.5*log(det(Sigma))+ ...
            + log(P_i(j));
    end
    Posterior = Posterior/max(Posterior);   %this will prevent numerical overflow
    
    Posterior = exp(Posterior)/sum(exp(Posterior));  %normalize
    [TempVal, TempInd] = max(Posterior);    %maximize over classes
    OptClass(i) = TempInd;
    MaxPosterior(i) = TempVal;
end

MP = max(MaxPosterior);
ind = find(MaxPosterior == MP);

% if length(ind) > 1
%     disp('there is a tie in the choice of subspace')
% end

indm = randperm(length(ind)); %break ties randomly

OptSubInd = ind(indm(1));
OptClass = OptClass(OptSubInd) - 1;


end

