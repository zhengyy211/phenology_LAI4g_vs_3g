function [data_lai, data_weights] = read_gimms_lai4g(dir_in, year)
% -------------------------------------------------------------------------
% Function:    read_gimms_lai4g
% Author:      Yaoyao Zheng
% Email:       zhengyaoyao@stu.pku.edu.cn
% MATLAB:      R2026a
%
% Description:
%   Reads and processes the GIMMS LAI4g semi-monthly dataset for a specified
%   year. The function:
%     1. Reads all 24 semi-monthly LAI files.
%     2. Converts fill values to NaN.
%     3. Applies the LAI scaling factor.
%     4. Extracts quality control (QC) flags.
%     5. Assigns quality weights based on the QC information.
%
% INPUTS:
%   dir_in        - Directory containing GIMMS LAI4g GeoTIFF files.
%   year          - Year of the data (e.g., 2001).
%
% OUTPUT:
%   data_lai      - Scaled LAI data (rows × columns × 24 semi-monthly periods).
%   data_weights  - Quality weights corresponding to data_lai.
%
% -------------------------------------------------------------------------


%% Step 1: Generate semi-monthly time sequence
months     = 1:12;
halves     = [1 2]; % 1 = first half-month, 2 = second half-month
[mm, hh]   = ndgrid(months, halves);
time_table = table(mm(:), hh(:), 'VariableNames', {'month', 'half_month'});
time_table = sortrows(time_table, {'month', 'half_month'});
n_time     = height(time_table);


%% Step 2: Read all LAI files
data_lai = cell(n_time, 1);
data_qc  = cell(n_time, 1);

for i_time = 1:n_time
    month = time_table.month(i_time);
    half = time_table.half_month(i_time);
    [data_lai{i_time}, data_qc{i_time}] = ...
        read_gimms_lai4g_hm(dir_in, year, month, half);
end

data_lai = cat(3, data_lai{:});
data_qc  = cat(3, data_qc{:});


%% Step 3: Convert QC flags to quality weights
% QC description:
% For 1982—2015, QC file is inherited from GIMMS NDVI3g product:
% 0: good quality for NDVI with no apparent problems;
% 1: NDVI from spline interpolation;
% 2: NDVI from seasonal profile due to possible snow or cloud cover;
% 3: missing value

% Default weight for low-confidence observations.
data_weights = ones(size(data_qc), 'single') * 0.2;

data_weights(data_qc == 0) =   1;  % Good-quality observation
data_weights(data_qc == 1) = 0.8;  % Spline interpolation
data_weights(data_qc == 2) = 0.5;  % Possible snow/cloud contamination
data_weights(data_qc == 3) = 0.2;  % Missing value

end


function [lai, qc] = read_gimms_lai4g_hm(dir_in, year, month, half_month)
% Reads one semi-monthly GIMMS LAI4g GeoTIFF file.
%
% INPUTS:
%  dir_in      - Input directory.
%  year        - Year of the data (e.g., 2001).
%  month       - Month (1-12).
%  half_month  - Half-month index (1 or 2).
%
% OUTPUT:
%   lai        - Scaled LAI values.
%   qc         - Quality control flags.

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


%% Step 3: Read LAI and QC layers
img = imread(file_path);
lai = single(img(:, :, 1));  % Extract LAI data
qc  = single(img(:, :, 2));  % Extract quality control flags


%% Step 4: Handle missing/invalid data
lai(lai == fill_val) = NaN;  % Replace fill value with NaN
qc(qc   == fill_val) = NaN;

lai = lai * lai_scale;  % Scale LAI values
end