function [lai, qc] = read_gimms_lai4g(dir_in, year, month, half_month)
% -------------------------------------------------------------------------
% Function:    read_gimms_lai4g
% Author:      Yaoyao Zheng
% Email:       zhengyaoyao@stu.pku.edu.cn
% Description: This function reads and processes the GIMMS LAI4g data for 
%              a specified year, month, and half-month period. It handles
%              missing or invalid data and scales the LAI values, and extracts the QC flags.
%
% INPUTS:
%   dir_in      - Directory containing the GIMMS LAI4g data files.
%   year        - Year of the data (e.g., 2001).
%   month       - Month of the data (1-12).
%   half_month  - Half-month period (1 or 2).
%
% OUTPUT:
%   lai         - Processed LAI data (scaled, missing values handled).
%   qc          - Quality control flags for the LAI data.
% 
% QC description:
%     0: good quality for NDVI with no apparent problems;
%     1: NDVI from spline interpolation;
%     2: NDVI from seasonal profile due to possible snow or cloud cover;
%     3: missing value

% -------------------------------------------------------------------------

%% Step 1: Define constants
fill_val  = 65535;  % Fill value for invalid data
lai_scale = 0.001;  % Scaling factor for LAI values


%% Step 2: Construct file name: 'GIMMS_LAI4g_V1.3_<year><month><half-month>.tif'
file_prefix = 'GIMMS_LAI4g_V1.3_';
file_name   = sprintf('%s%d%02d%02d.tif', file_prefix, year, month, half_month);
file_path   = fullfile(dir_in, file_name);

if ~isfile(file_path)
    error('File %s does not exist.', file_path);
end


%% Step 3: Read data
img = imread(file_path);
lai = single(img(:, :, 1));  % Extract LAI data
qc  = single(img(:, :, 2));  % Extract quality control flags


%% Step 4: Handle missing/invalid data
lai(lai == fill_val) = NaN;  % Replace fill value with NaN
lai = lai * lai_scale;       % Scale LAI values
qc(qc == fill_val) = NaN;    % Replace fill value in QC flags with NaN

end