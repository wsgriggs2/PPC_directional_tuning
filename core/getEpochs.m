function [fix, mem, mov] = getEpochs(behavior)
% [fix, mem, mov] = getEpochs(behavior)
% returns the epoch/phase timing of the trial based on the behavior
% example:
% behavior is a 1xnTrials struct with fields, e.g. trialID, trialstart,...
% fix = [-5.2 0]
% mem = [0  4.9]
% mov = [4.9 12]
% note: all values are in seconds relative to initial/memory cue onset


% only look at successful trials
success = cell2mat({behavior.success});
behavior(~success) = [];

% get vector of times for each phase
trialstart = parseTime(behavior, 'trialstart');
fixationHold = parseTime(behavior, 'fixationhold');
cue = parseTime(behavior, 'cue');
try
    memoryStart = parseTime(behavior, 'memory');
    AlignToTime = 'memory';
catch
   warning('No memory period. Aligning to cue period instead.');
   AlignToTime = 'cue';
end
        
goCue = parseTime(behavior,'target_acquire');
reward = parseTime(behavior, 'reward');
iti = parseTime(behavior, 'iti');

% realign to memory or cue start
switch AlignToTime
    case 'memory'
        fixationHold = fixationHold - memoryStart;
        cue = cue - memoryStart;
        goCue = goCue - memoryStart;
        reward = reward - memoryStart;
    case 'cue'
        fixationHold = fixationHold - cue;
        goCue = goCue - cue;
        reward = reward - cue;
end

% set periods
fix = [median(fixationHold) 0];
mem = [0 median(goCue)];
mov = [median(goCue) max(reward)];

iti_length = mean(trialstart(2:end)-iti(1:end-1), 'omitnan');

end