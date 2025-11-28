function [matrix_lai, mask_spike_all] = remove_spike(matrix_lai, matrix_qc, good_qc, matrix_thermal_gs)
% -------------------------------------------------------------------------
% Function:    remove_spike
% Author:      Yaoyao Zheng
% Email:       zhengyaoyao@stu.pku.edu.cn
% Description: This function removes spikes from LAI data by comparing 
%              each observation with a three-point moving average. 
%
%              **Spike correction is applied only during the non-growing season** 
%              (i.e., when temperature is below 0°C). During the thermal growing 
%              season (T > 0°C), no correction is applied.
%
%              **Spike Definition**:
%              A spike is defined as a data point where the absolute difference 
%              between the observed LAI value and the three-point moving average 
%              exceeds 25% of the maximum value of that pixel across the entire 
%              time series.
%   
% INPUTS:
%   matrix_lai        - Matrix containing the LAI data [pixels x time].
%   matrix_qc         - Matrix containing the quality control flags [pixels x time].
%   good_qc           - The quality control values that represent "good" data.
%   matrix_thermal_gs - Logical matrix for the thermal growing season [pixels x time].
%
% OUTPUTS:
%   matrix_lai        - Matrix containing the corrected LAI data.
%   mask_spike_all    - Logical matrix indicating the locations of detected spikes.
% -------------------------------------------------------------------------

% Compute pixel-level baseline
% (minimum and maximum of good-quality observations)
data_good = matrix_lai;
data_good(~ismember(matrix_qc, good_qc)) = NaN;
data_min = min(data_good, [], 2, 'omitmissing');
data_max = max(data_good, [], 2, 'omitmissing');
clearvars data_good
data_min = max(data_min, 0);  % Ensure the minimum value is non-negative
data_min = repmat(data_min, 1, size(matrix_lai, 2));

done     = false;
iter     = 0;
max_iter = 3;   % Maximum number of iterations for spike removal
mask_spike_all = false(size(matrix_lai));  % Initialize mask for spikes

while ~done && iter < max_iter
    iter = iter + 1;

    % Compute the three-point moving average
    lai_mavg = movmean(matrix_lai, 3, 2, 'omitmissing');

    % Detect spikes: values that deviate more than 25% of the maximum value
    mask_spike = abs(matrix_lai - lai_mavg) > 0.25 .* data_max;
    
    % Apply correction only during thermal non-growing season (T < 0°C)
    mask_spike = mask_spike & (matrix_thermal_gs == false);

    if ~any(mask_spike, 'all')
        done = true;
    else
        % Correct the detected spikes by setting them to the minimum of good-quality data
        matrix_lai(mask_spike) = data_min(mask_spike);
        % Update the global spike mask
        mask_spike_all = mask_spike_all | mask_spike;
    end
end

end