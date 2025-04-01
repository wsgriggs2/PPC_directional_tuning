function [percent_depth, actual_depth] = get_percent_depth(region_mask)
    min_depth = find(any(region_mask, 2), 1);
    max_depth = find(any(region_mask, 2), 1, 'last');
    
    % Find y-coordinate of each voxel in the mask
    [y_coords, ~] = find(region_mask);
    
    % Convert y-coordinate into percent of total LIP depth
    percent_depth = (y_coords - min_depth)/(max_depth - min_depth);
    actual_depth = y_coords - min_depth;
end