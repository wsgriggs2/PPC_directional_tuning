function behavior_out = adjust_behavior_time(behavior_in, statenames, time_amount)
% Function to add `time` to every event time in behavior_in

behavior_out = behavior_in;

for field_ind = 1:length(statenames)
    columnTitle = statenames{field_ind};
    

    
    %% Extract time in seconds and subtract correction
    % Extract time

    event_time_orig = parseTime(behavior_in, columnTitle);

   
   % Remove NaN values
   event_time_orig(isnan(event_time_orig)) = [];
   event_time_corrected = event_time_orig + time_amount;
   
   %% Rebuild the strings - opposite of parseTime
   % Find which cells are empty in the struct
   temp_cell = {behavior_in.(columnTitle)}';
   empty_ind = cellfun(@isempty, temp_cell);
   
   % Pre-allocate space
   numTimeArray = zeros(size(event_time_corrected, 1), 4);
   
   % Would be simpler to use `seconds` function, but for some reason, it
   % will not work when it is in this function...
   negative_ind = event_time_corrected < 0;
   
   %Convert back into hours, minutes, seconds, and milliseconds
   event_time_corrected = abs(event_time_corrected);
   numTimeArray(:, 1) = floor(event_time_corrected/3600);
   numTimeArray(:, 2) = floor((event_time_corrected - numTimeArray(:, 1)*3600)/60);
   numTimeArray(:, 3) = floor(event_time_corrected - numTimeArray(:, 1)*3600 - numTimeArray(:, 2)*60);
   numTimeArray(:, 4) = rem(event_time_corrected, 1);
   

   if max(numTimeArray(:, 1)) >= 100
       % String format will be HHH:MM:SS,FFF
       hours = num2str(numTimeArray(:, 1), '%03d');
   else
       % String format will be HH:MM:SS,FFF
       hours = num2str(numTimeArray(:, 1), '%02d');
   end
   minutes = num2str(numTimeArray(:, 2), '%02d');
   seconds = num2str(numTimeArray(:, 3), '%02d');
   milliseconds = num2str(numTimeArray(:, 4), '%0.3f');
   
   colons = repmat(char(':'), size(numTimeArray, 1), 1);
   commas = repmat(char(','), size(numTimeArray, 1), 1);
   
   charTime = [hours colons minutes colons seconds commas milliseconds(:, 3:5)];
   
   i = 0;
   for ind = find(~empty_ind)'
    i = i+1;
    if negative_ind(i)
        behavior_out(ind).(columnTitle) = ['-' charTime(i, :)];
    else
        behavior_out(ind).(columnTitle) = charTime(i, :);
    end
   end
    
end