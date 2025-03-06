function [filt_data, filt_time] = std_filt(data, time, std_fact)

% modifications made from lab 1 are to accept data in a sequential column
% format instead of data being encoded in sequential rows and to simplify
% clamp the values if necessary (as data is needed continuously at the
% original timestep to ensure ekf works properly)

[raw_var, raw_mean] = var(data,0,2); % sample statistics
raw_std = sqrt(raw_var);

bounds_offset = std_fact * raw_std;

% removes any data outside of the bounds formed by standard deviation
outliers = abs(data - raw_mean) > bounds_offset;
valid_points = all(~outliers,1); % ensures all rows have to be valid

filt_data = data(valid_points);
filt_time = time(valid_points);

end