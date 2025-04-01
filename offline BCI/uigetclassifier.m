function [classifier_string, validation_string] = uigetclassifier()
% this function asks the user what type of classifier and what type of
% cross validation they would like to use and returns strings for each.
% Those strings are usually fed to classifyDoppler.m
% Author: Sumner Norman
% Date: 2019
validationTypes = {'kFold','leaveOneOut'};
classifierTypes = {'CPCA+LDA', 'PCA+LDA'};

[classifierIdx,~] = listdlg('PromptString','Select Classifier Type:',...
    'SelectionMode','single',...
    'ListString',classifierTypes);

[validationIdx,~] = listdlg('PromptString','Select Validation Type:',...
    'SelectionMode','single',...
    'ListString',validationTypes);

if isempty(classifierIdx*validationIdx)
    warning('invalid selections'), return,
end

classifier_string = classifierTypes{classifierIdx};
validation_string = validationTypes{validationIdx};
end