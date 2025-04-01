function [pred, model] = classifyDoppler(trainData, trainLabels, testData, varargin)

% class = classifyDoppler(trainData, trainLabels, method, testData)
% classifies doppler data using the training data and labels provided
%
% varargin:
%
% trainData --  training data, size nImagesL&R x (nPixels*yPixels)
% trainLabels -- training labels, size nImagesL&R x 1
% testData: testingData of size nImagesToTest x (nPixels*yPixels)
% method -- string of method type. Accepted values are:
%   'CPCA+LDA'
%   'PCA+LDA'
% subspace -- optional - you can pass in the subspace you want it to use
%             to classify. Note that this only works with CPCA+LDA!
%             Typically, you would pass in DRmatC which is a cell array of
%             2D doubles.
% variance_to_keep -- optional - for the PCA dimensionality reduction
%                     methods. Not currently used for cPCA.
%
% varargin:
% 'cPCA_m' - Scalar positive value; determines dimension of each cPCA
% subspace. Max is (# of classes - 1). Default is 1.
%
% varargout:
%
% pred: a vector value of predicted class for each replicate (previously
% class but changed to be a more general variable)
% model:  model used to generate predictions (previously subspace)
%
% Note that the subspace varargin/varargout are really only applicable for
% the CPCA+LDA cases. This is probably bad practice since we have arguments
% that don't make sense for all of the classification types here. We should
% probably pass out a generic structure that works for all classifier types

%% Variable argument input parser
p = inputParser;
p.addParameter('m', 1, @isscalar);
p.addParameter('method','CPCA+LDA');
p.addParameter('model', NaN);
p.addParameter('variance_to_keep', 95);

p.parse(varargin{:});
VariableInputs = p.Results;
m = VariableInputs.m;
method = VariableInputs.method;
model = VariableInputs.model;

% remove nans if appropriate
nan_i = isnan(trainLabels);
trainLabels(nan_i,:) = [];
trainData(nan_i,:) = [];

%--------------------------------------------- PCA+LDA ---------------------------------------------
if strcmp(method,'PCA+LDA')
    
    if ~isstruct(model)
        % Needs to be nested, because cannot use isnan on a struct.
        if isnan(model)
            % PCA
            
            [pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(trainData);
            explainedVarianceToKeepAsFraction = VariableInputs.variance_to_keep/100;
            if VariableInputs.variance_to_keep < 100
                numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
                pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
                trainPredictors = pcaScores(:,1:numComponentsToKeep);   % (PCA transformed trainData)
                
            else
                trainPredictors = pcaScores;   % (PCA transformed trainData)
            end
            
            
            MdlLinear = fitcdiscr(trainPredictors,trainLabels);
            
            model = struct();
            model.pcaCenters = pcaCenters;
            model.pcaCoefficients = pcaCoefficients;
            model.MdlLinear = MdlLinear;
        end
    else
        pcaCenters = model.pcaCenters;
        pcaCoefficients = model.pcaCoefficients;
        MdlLinear = model.MdlLinear;
        
        
    end
    testDataPCA = (testData - pcaCenters) * pcaCoefficients;

    pred = predict(MdlLinear, testDataPCA);


    %-------------------------------------------- CPCA+LDA ---------------------------------------------
elseif strcmp(method,'CPCA+LDA')
    
    if ~iscell(model) && isnan(model)
        % create transformations
        DRmatC = trainCPCA(trainData, trainLabels, 'm', m);
    else
        DRmatC = model;
    end
    
    % make sure it didn't pick an empty subspace (happening
    % sometimes as of 9/6/2018 - this is a workaround)
    emptyInds = cellfun('isempty',DRmatC);
    if any(emptyInds)
        warning('There are empty subspaces here!!!')
        firstNonEmpty = find(~emptyInds); firstNonEmpty = firstNonEmpty(1);
        DRmatC{emptyInds} = DRmatC{firstNonEmpty};
    end
    
    % determine the optimal subspaces (for each test trial)
    nTrialsTest = size(testData,1); % number
    S = zeros(1,nTrialsTest);       % which subspace was best
    for trial = 1:nTrialsTest
        [S(trial),~] = choose_subspace(testData(trial,:),...
            trainData,trainLabels,DRmatC,'unbiased');
    end
    
    % classify for all trials in each subspace:
    pred = zeros(nTrialsTest, 1);
    for s = unique(S)
        ind = find(S == s);
        testFeatures = testData(ind,:) * DRmatC{s};
        trainFeatures = trainData * DRmatC{s};
        pred(ind) = classify(testFeatures, trainFeatures, trainLabels,...
            'linear','empirical'); %#ok<*AGROW>
    end
    
    % rename DRmatC for export
    model = DRmatC;
else
    error('This decoding method is not currently supported!');
end