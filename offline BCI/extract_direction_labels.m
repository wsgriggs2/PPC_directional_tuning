function [TargetPosInd, train_labels, label_strings] = extract_direction_labels(behavior)
target = cell2mat({behavior.target});
targetPos = vertcat(behavior.targetPos);
if iscell(targetPos)
    targetPos = cell2mat(targetPos);
end
[angle,distance] = cart2pol(targetPos(:,1),targetPos(:,2)); % convert to polar coordinates
angle(angle<0) = angle(angle<0)+2*pi; %Convert to have only positive angles.

UniqueAngles = unique(angle);
TargetPosInd = zeros(length(target),1);
for position = 1:length(UniqueAngles)
    TargetPosInd(angle==UniqueAngles(position)) = position;
end

% Decide if we want to use the center positions.
useMiddlePosition = true;

% create the training data labels (2 x n_trials) - x, y columns
% For horizontal and vertical axes, use the following labels
% 1 = Negative (down or left)
% 2 = At center (along vertical or horizontal axis)
% 3 (or 2 if not using middle points) = Positive (up or right)
if useMiddlePosition
    label_strings = {'Negative','Center','Positive'};
else
    label_strings = {'Negative','Positive'};
end


train_labels = NaN(size(targetPos,1),2);

for dimension = 1:2
    train_labels(targetPos(:, dimension) < 0, dimension) = 1;
    if useMiddlePosition
        train_labels(targetPos(:, dimension) == 0, dimension) = 2;
        train_labels(targetPos(:, dimension) > 0, dimension) = 3;
    else
        % As of 1/6/21, we are ignoring these middle points. Leaving this
        % code in case it is useful. It randomly assigns the middle points
        % to either class 1 or 2.
        % train_labels(targetPos(:, dimension) == 0, dimension) = randi(2, nnz(targetPos(:, dimension) == 0),1);
        
        % If ignoring midle points, then there are only two classes
        train_labels(targetPos(:, dimension) > 0, dimension) = 2;
    end
    

end