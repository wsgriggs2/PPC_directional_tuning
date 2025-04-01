function angiogram = makeAngiogram(DopplerImage, varargin)

%% Variable inputs
p = inputParser;
p.addOptional('root',3)
p.parse(varargin{:});
inputs = p.Results;
%% Make the angiogram

if ndims(DopplerImage) == 4
    angiogram = nthroot(squeeze(mean(mean(DopplerImage,3),4)),inputs.root);

elseif ndims(DopplerImage) == 3
    angiogram = nthroot(squeeze(mean(DopplerImage, 3)), inputs.root);
elseif ndims(DopplerImage) == 2
    angiogram = nthroot(DopplerImage, inputs.root);
end
lowerCutoff = quantile(angiogram(:),0.01);
angiogram(angiogram<lowerCutoff) = lowerCutoff;
    
end