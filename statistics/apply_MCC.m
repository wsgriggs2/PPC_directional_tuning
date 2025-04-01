function pvalues_corrected = apply_MCC(pvalues, varargin)
%% This function applies multiple comparison correction
% Inputs:
%   pvalues - Any size array
%   'method' - 'FDR' or 'Bonferroni'; Would be nice to add FWE in the
%               future
% Outputs:
%   pvalues_corrected - Corrected p-values in same shape as input array
%
% Notes:
% 1. Currently will not handle NaN values. This is by design. No NaN
%    values should ever be passed to this function.
% 
% Written February 2, 2022 Whitney Griggs
%% Parse inputs
p = inputParser;
addOptional(p, 'method',  'FDR');
addOptional(p, 'verbose', true);
parse(p, varargin{:});
result = p.Results;

%% Check for NaN values
if any(isnan(pvalues), 'all')
   error('This function will not handle NaN values correctly. Please pass pvalues within [0 1] to this function'); 
end

%% Apply multiple comparison correction method of choice
switch result.method
    case 'FDR'
        % Reshape into 1D array of values
        ndims = length(size(pvalues));
        resize_array = ones(1, ndims);
        resize_array(1) = numel(pvalues);
        p_1d = reshape(pvalues, resize_array);
        
        % Apply FDR correction
        [pvalues_corrected_1d] = mafdr(p_1d, 'BHFDR', true, 'Showplot', result.verbose);
        
        % Reshape into original size
        pvalues_corrected = reshape(pvalues_corrected_1d, size(pvalues));

    case 'bonferroni'
        % Apply bonferroni correction
        pvalues_corrected = pvalues./numel(pvalues);
    otherwise
        error('"%s" is not a currently supported method. Please choose from "FDR" or "bonferroni"', result.method);
end
