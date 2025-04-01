function numTime = parseTime(behavior_struct, columnTitle) 
%% Converts the character arrays in the behavior data to numbers
% Moved from being a subfunction in parsefUSBehavior to being an independent
% function due a restructuring of parsefUSBehavior requiring this function in
% multiple different functions
% 
% Should replace `parseTime` that was previously in use (core/parseTime.m).
%
% Written by Sumner
% Modified by Whitney to make more general and handle empty time cells


    temp_cell = {behavior_struct.(columnTitle)}';
    empty_ind = cellfun(@isempty, temp_cell);
    negative_ind = false(size(empty_ind));
    negative_ind(~empty_ind) = contains(temp_cell(~empty_ind), '-');
    
    % The first entry should always be the largest since it will the the
    % length of the standard character + a minus sign if it exists. This may break if
    % the hours ever gets to be 100 or greater though. 
    number_chars = size(temp_cell{1}, 2);
    if any(negative_ind)
        chartime(~empty_ind & ~negative_ind, 2:number_chars) = cell2mat(temp_cell(~empty_ind & ~negative_ind));
        chartime(~empty_ind & negative_ind, :) = cell2mat(temp_cell(~empty_ind & negative_ind));
    else
        chartime = cell2mat(temp_cell(~empty_ind));
    end

    % extract hours/min/sec/ms from chartimes
    %#ok<*ST2NM>
    
    % For flexibility, we are parsing this a single time string at a time.
    % Previously, we could use the assumption of HH:MM:SS,FFF to vectorize
    % the code, but with the advance to continuous data, there is (remote)
    % possibility of HHH:MM:SS,FFF.
    
    numTimeArray = zeros(size(chartime, 1), 3);
    
    for i = 1:size(chartime, 1)
        temp_chartime = string(chartime(i, :));
        colon_ind = strfind(temp_chartime, ':');
        comma_ind = strfind(temp_chartime, ',');
        negative_ind(i) = contains(temp_chartime, '-');
        numTimeArray(i, 1) = str2num(chartime(i,1:colon_ind(1)-1));   % Hours
        numTimeArray(i, 2) = str2num(chartime(i,colon_ind(1)+1:colon_ind(2)-1));   % Minutes
        numTimeArray(i, 3) = str2num(chartime(i,colon_ind(2)+1:comma_ind-1));   % seconds
        numTimeArray(i, 3) = numTimeArray(i, 3) + 1e-3*str2num(chartime(i,comma_ind+1:comma_ind+3)); % milliseconds
    end


    % Keep same size as input
    % Using max(size) because this seems a little more general than
    % assuming behavior struct is always 1 x N.
    numTime = NaN(max(size(behavior_struct)), 1);
    
    % convert time to [s]
    numTime(~empty_ind) =   numTimeArray(:,1)*3600 + numTimeArray(:,2)*60+ numTimeArray(:,3);
    numTime(negative_ind) = numTime(negative_ind)*-1;
end