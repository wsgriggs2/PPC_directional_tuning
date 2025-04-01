function DRmatC = trainCPCA(trainData, trainLabels, varargin)
% trainCPCA will find the transformation subspaces (one for each class)
% using classwise principal component analysis (CPCA). There are several
% options defining specifics of the CPCA algorithm. One of the more 
% important is the DFEmethod, which specifies the discriminant feature 
% extraction method to use. For example, this can be toggled between LDA 
% and IDA, etc.
%
% CPCA Options:
% m (1 x 1) is the final dimension of the extracted subspaces
%
% varargin:
% m - Scalar positive value; determines dimension of subspaces. Max is (# of
% classes - 1). Default is 1.
%
% Prior (string) specifies the way prior (class) probabilities are
% calculated in the method. 'empirical' means that the class relative 
% frequencies are preserved; in our example [0.5 0.5], since there are 20 
% instances of dopL and 20 instances of dopR.
%
% EvalKeep (cell) specifies which eigenvalues/eigenvectors to keep.
%
% DFEmethod (string) specify discriminant feature extraction method to use


%% Variable argument input parser
p = inputParser;
p.addParameter('m', 1, @isscalar)

p.parse(varargin{:});
VariableInputs = p.Results;
m = VariableInputs.m;

%% CPCA setup
prior = 'uniform';
EvalKeep = {'mean'};
DFEmethod = 'lda';

%% remove nans (redundant/just in case - usually done by classifyDoppler)
nan_i = isnan(trainLabels); 
trainLabels(nan_i,:) = [];
trainData(nan_i,:) = [];

%% we can now run CPCA 
% DRmatC = dataproc_func_cpca(trainData,trainLabels,m,Prior,EvalKeep,DFEmethod)
%
% DRmatC (1 x C) cell, whrere C is the number of classes.
%
% m (1 x 1) is the final dimension of the extracted subspaces
%
% Prior (string) specifies the way prior (class) probabilities are
% calculated in the method. 'empirical' means that the class relative 
% frequencies are preserved; in our example [0.5 0.5], since there are 20 
% instances of dopL and 20 instances of dopR.
%
% EvalKeep (cell) specifies which eigenvalues/eigenvectors to keep.
%
% DFEmethod (string) specify discriminant feature extraction method to use

DRmatC = dataproc_func_cpca(trainData, trainLabels, m, prior, EvalKeep, DFEmethod);

% The function returns two 1D subspaces: DRmatC{1} and DRmatC{2}. These 
% subspaces (one for each class) appear due to nonlinearity (piecewise 
% linearity) of CPCA. Each subspace can be seen as feature extraction 
% mapping from the original N-dimensional data space into 1D feature 
% subspace. They also have a physical interpretation as they point to the 
% areas of image (pixels) that encode the differences between dopL and dopR. 

end

