function dopOut = preProcess(dopIn, varargin)
% dopOut = preProcess(dopIn, vargin)
%
% preProcess will perform basic pre-processing steps on the raw doppler
% data of the form: [yPixels x xPixels x timeWindows x trials]
% Good practice would require that all doppler data is passed together,
% regardless of target. For example, pass iDop, but not dopIn.
%
% important: preProcess is usually called immediately after loading a
% session's data using loadSessionData
%
% optional vargin:
% 'timeGain',               bool
% 'zscore',                 bool
% 'spatialFilter',          cell array; Parameters that will be eventually
%                           passed to FSPECIAL.
% 'verbose',                bool
%
% vargout
% dopOut is of the same form as dopIn, but with preprocessing steps

%% handling varargin
p = inputParser;
p.addOptional('timeGain',false,@islogical);
p.addOptional('zScore',false,@islogical);
p.addOptional('spatialFilter',{'disk', 2});
p.addOptional('verbose',false,@islogical);
p.parse(varargin{:});
sets = p.Results;
if sets.verbose
    disp(sets)
end

%% main function is here (function calls)
if sets.timeGain, dopIn = timeGainCompensation(dopIn,sets); end
if sets.zScore, dopIn = dopZ(dopIn); end
if any(~cellfun(@isempty, sets.spatialFilter)), dopIn = spatialFilter(dopIn, sets.spatialFilter); end

dopOut = dopIn;

%% spatial filter
    function dataOut = spatialFilter(dataIn, filter_options)
        % Apply specified spatial filter
        % First entry in `filter_options` cell defines the spatial filter
        % Subsequent entries define the parameters for that specified
        % filter.
        %
        % For 'disk' - specify filter size
        % For 'gaussian - specify filter size and filter sigma
        [~, ~, n_windows, n_trials] = size(dataIn);
        switch filter_options{1}
            case 'disk'
                disk_size = filter_options{2};
                if disk_size
                    if disk_size && islogical(disk_size)
                        % For backwards compatibility. True used to always
                        % mean `disk_size=2`.
                        disk_size = 2;
                    end
                    h = fspecial(filter_options{1}, disk_size); % disk size defined here
                    dataOut = filter_doppler_2D(dataIn, h, n_windows, n_trials);
                else
                    dataOut = dataIn;
                end
            case 'gaussian'
                h = fspecial(filter_options{1}, filter_options{2}, filter_options{3}); % disk size defined here
                dataOut = filter_doppler_2D(dataIn, h, n_windows, n_trials);
            otherwise
                error('This filter type has not been implemented');
        end



    end

%% time gain compensation
    function dataOut = timeGainCompensation(dataIn, res)
        % getting useful dimensions
        [n_depth, n_width, n_windows, n_trials] = size(dataIn);
        depthInds = (1:n_depth)';
        indicesForModeling = round(0.1*n_depth):n_depth; % discards supra-dural info

        % get mean signal (empirical)
        dataIn_flat = reshape(dataIn, n_depth, n_width, []);
        depthMeanSignal = squeeze(mean(mean(dataIn_flat,2),3)); % 1D depth array (actual)

        % create the exponential model for depth attenuation/time-gain-comp
        f = fit(depthInds(indicesForModeling) , double(depthMeanSignal(indicesForModeling)), 'exp1');
        depthMean = f.a*exp(f.b*depthInds);

        % plot the result
        if res.verbose
            figure();
            plot(depthInds, depthMean, 'k', 'lineWidth', 2);
            hold on; plot(depthInds, depthMeanSignal, '.r')
            legend('depth attenuation model','actual mean activation')
            xlabel('depth (pixels)'); ylabel('raw signal [a.u.]')
        end

        % normalization happens here
        dataOut = dataIn./repmat(depthMean, 1, n_width, n_windows, n_trials);
        dataOut = dataOut*100;  % converts to percentage of mean signal
    end

%% z score each voxel
    function dataOut = dopZ(dataIn)
        [n_depth, n_width, n_windows, n_trials] = size(dataIn);

        % flatten the data into a voxel-column for zscoring
        dataOut = permute(dataIn, [2, 1, 4, 3]);% [n_width, n_depth, n_trials, n_windows]);
        dataOut = reshape(dataOut, [n_depth*n_width, n_windows*n_trials]);
        % perform the zscore
        dataOut = zscore(dataOut')';
        % reshape back to original data dims
        dataOut = reshape(dataOut, [n_width, n_depth, n_trials, n_windows]);
        dataOut = permute(dataOut, [2, 1, 4, 3]);
    end

%% Helper functions
    function dataOut = filter_doppler_2D(dataIn, h, n_windows, n_trials)
        dataOut = NaN(size(dataIn));
        for window = 1:n_windows
            for trial = 1:n_trials
                dataOut(:, :, window, trial) = ...
                    filter2(h, squeeze(dataIn(:, :, window, trial)));
            end
        end

    end
end