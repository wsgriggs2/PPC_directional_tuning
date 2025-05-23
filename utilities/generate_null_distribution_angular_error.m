function [null_angular_error, null_mean_angular_error] = generate_null_distribution_angular_error(nClasses, num_repetitions)
% generate_null_distribution_angular_error  Generate the null distribution
% for angular error by randomly sampling from a uniform distribution of the
% possible directions.
%
% INPUTS:
%   nClasses:                   double; 2 or 8
%   num_repetitions:            scalar double; (n); How many replicates?
%      
% OUTPUTS:
%   null_angular_error:         (n x m) double; Randomly sampled angular
%                               error drawn from uniform distribution for
%                               each trial. n = num_replicates. m =
%                               num_trials
%   null_mean_angular_error:    (n x m) double; Trial-wise iterative mean
%                               of the randomly angular error. 
%                               n = num_replicates. m = num_trials

if ~exist('nClasses', 'var')
    nClasses = [8];
end

if ~exist('num_repetitions', 'var')
    num_repetitions = 1000;
end

for i = 1:length(nClasses)
    switch nClasses(i)
        case 8
            % Map from predicted number to direction
            angle_lookup_table = 0:pi/4:7*pi/4;
            angle_lookup_table = angle_lookup_table([6 7 8 5 1 4 3 2]);
        case 2
            % Map from predicted number to direction
            angle_lookup_table = [0 pi];
        otherwise
            error('%d number of classes is not currently supported.', nClasses);
    end
    
    % Very generous upper bound of number of trials ever expected in single session
    nTrials = 2000; 
    
    % Pre-allocate space
    null_mean_angular_error = NaN(num_repetitions, nTrials);
    
    
    % Calculate angular error for each trial
    % Randomly sample form uniform distribution
    random_direction = randi(nClasses(i), num_repetitions, nTrials);
    null_angular_error = abs(angdiff(angle_lookup_table(random_direction), zeros(num_repetitions, nTrials)));
    
    ppm = ProgressBar(nTrials*num_repetitions, taskname='Generating null distribution');
    parfor trial = 1:nTrials
        % Calculate iterative mean as we add trials.
        for rep = 1:num_repetitions
            null_mean_angular_error(rep, trial) = mean(null_angular_error(rep, 1:trial));
            ppm.count();
        end
    end
    delete(ppm);
    
    %% Save the null distribution
    save(sprintf('angular_error_null_distribution_%dtgts_%dreps.mat', nClasses(i), num_repetitions), 'null_angular_error', 'null_mean_angular_error');
    
end