function circular_stats = get_circular_statistics(theta, rho, varargin)

%% Compute multiple circular statistic summary values
% Inputs:
%   theta - m x n; where m is num_of_pixels and n is the number of
%           different directions. Direction of each target in radians;
%           Technically, only first row is needed
%   rho -   m x n; where m is number of pixels and n is number of different
%           directions. Each row represents the response of that voxel to
%           each of the n directions.
%
% Variable inputs
%   `bin_spacing` - Default is pi/4 (45d); needs to be equal bin spacing.
%   'x_pixels' - Default is 128; Assumes 128 element transducer
%
% Output
%   circular_stats - struct with the following summary statistics generated
%                    for each voxel - variance, std, skewness, kurtosis,
%                    parameter fit for Von Mise function. Some of these
%                    have multiple possible versions (indicated by `_alt`
%                    suffix.
%
% Written by Whitney Griggs - January 2023.

%% Input parser
p = inputParser;
addOptional(p, 'bin_spacing',  pi/4); % In radians
addOptional(p, 'x_pixels', 128); %Assumes 128 element transducer
addOptional(p, 'verbose', false);
parse(p, varargin{:});
inputs = p.Results;

bin_spacing = inputs.bin_spacing;
x_pixels = inputs.x_pixels;

%% Get circular statistics for each voxel
num_pixels = size(rho, 1);

% What are the positions of the targets?
alpha = theta(1, :);

[r, thetahat, kappa, S, s, cstd, cstd0, b, b0,...
    k, k0, pk, loc, width, prom, n_pks] = deal(NaN(num_pixels, 1));

if inputs.verbose
    wb = waitbar(0, 'computing circular statistics');
end
for i = 1:num_pixels
    response = rho(i, :);
    
    % If any response is less than zero, shift so that number is the
    % minimum.
    if any(response<0)
        response_offset = response - min(response);
    else
        response_offset = response;
    end
    r(i) = circ_r(alpha', response_offset', bin_spacing);
    [thetahat(i), kappa(i)] = circ_vmpar(alpha', response_offset', bin_spacing);
    [S(i), s(i)] = circ_var(alpha', response_offset', bin_spacing);
    [cstd(i), cstd0(i)] = circ_std(alpha', response_offset', bin_spacing);
    [b(i), b0(i)] = circ_skewness(alpha', response_offset');
    [k(i), k0(i)] = circ_kurtosis(alpha', response_offset');
    
    
    % For findpeaks, won't find peaks at edges, so need to add one entry to
    % allow for peaks at edges
    [pks, locs, widths, proms] = findpeaks(response_offset([end 1:end 1]), ...
        [alpha(1)-diff(alpha(1:2)) alpha alpha(end)+diff(alpha(end-1:end))], ...
        'WidthReference', 'halfheight', ...
        'SortStr', 'descend', ...
        'Annotate', 'extents');
    n_pks(i) = length(pks);
    if isempty(pks)
        if inputs.verbose
            fprintf('No peaks found\n');
        end
    else
        pk(i) = pks(1);
        loc(i) = locs(1);
        width(i) = widths(1);
        prom(i) = proms(1);

        

        % Handle edge case that loc is negative or greater than end alpha value
        if loc(i)<alpha(1)
            loc(i) = alpha(end);
        elseif loc(i)>alpha(end)
            loc(i) = alpha(1);
        end
    end

    if inputs.verbose
        waitbar(i/num_pixels);
    end
end

% Convert all locations to be within [-pi pi]
greaterThanPIind = loc > pi;
lessThanNegPIind = loc < - pi;

loc(greaterThanPIind) = loc(greaterThanPIind) - 2*pi;
loc(lessThanNegPIind) = loc(lessThanNegPIind) + 2*pi;

if inputs.verbose
    close(wb);
end

%% Transform each variable back into map
r_map = reshape(r, [], x_pixels);
thetahat_map = reshape(thetahat, [], x_pixels);
kappa_map = reshape(kappa, [], x_pixels);
S_map = reshape(S, [], x_pixels);
s_map = reshape(s, [], x_pixels);
cstd_map = reshape(cstd, [], x_pixels);
cstd0_map = reshape(cstd0, [], x_pixels);
b_map = reshape(b, [], x_pixels);
b0_map = reshape(b0, [], x_pixels);
k_map = reshape(k, [], x_pixels);
k0_map = reshape(k0, [], x_pixels);
pks_map = reshape(pk, [], x_pixels);
locs_map = reshape(loc, [], x_pixels);
width_map = reshape(width, [], x_pixels);
prom_map = reshape(prom, [], x_pixels);
n_pks_map = reshape(n_pks, [], x_pixels);

%% Put all into single struct
circular_stats.r = r_map;
circular_stats.thetahat = thetahat_map;
circular_stats.kappa = kappa_map;
circular_stats.var = S_map;
circular_stats.var_alt = s_map;
circular_stats.circ_std = cstd_map;
circular_stats.circ_std_alt = cstd0_map;
circular_stats.skewness = b_map;
circular_stats.skewness_alt = b0_map;
circular_stats.kurtosis = k_map;
circular_stats.kurtosis_alt = k0_map;
circular_stats.pks = pks_map;
circular_stats.locs = locs_map;
circular_stats.width = width_map;
circular_stats.prom = prom_map;
circular_stats.n_pks = n_pks_map;