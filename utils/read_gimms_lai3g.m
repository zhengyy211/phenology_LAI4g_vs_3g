function [data_lai, data_weights] = read_gimms_lai3g(dir_in, year)
% -------------------------------------------------------------------------
% Function:    read_gimms_lai3g
% Author:      Yaoyao Zheng
% Email:       zhengyaoyao@stu.pku.edu.cn
% MATLAB:      R2026a
%
% Description:
%   Reads and processes the GIMMS LAI3g semi-monthly dataset for a specified
%   year. The function:
%     1. Reads all 24 semi-monthly LAI files.
%     2. Converts fill values to NaN.
%     3. Applies the LAI scaling factor.
%     4. Extracts quality control (QC) flags.
%     5. Assigns quality weights based on the QC information.
%
% INPUTS:
%   dir_in        - Directory containing GIMMS LAI3g GeoTIFF files.
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
        read_gimms_lai3g_hm(dir_in, year, month, half);
end

data_lai = cat(3, data_lai{:});
data_qc  = cat(3, data_qc{:});


%% Step 3: Convert QC flags to quality weights
% QC description:
% Quality flag was embedded in the data. 
% Quality flag can be calculated by: flag = data - floor(data/10)*10
% flag == 0: No apparent issues (good value)
% flag == 1: Retrieved from spline interpolation
% flag == 2: Retrieved from seasonal profile (possible snow/cloud)

% Default weight for low-confidence observations.
data_weights = ones(size(data_qc), 'single') * 0.2;

data_weights(data_qc == 0) =   1;  % Good-quality observation
data_weights(data_qc == 1) = 0.8;  % Spline interpolation
data_weights(data_qc == 2) = 0.5;  % Possible snow/cloud contamination

end


function [lai, qc] = read_gimms_lai3g_hm(dir_in, year, month, half_month)
% Reads one semi-monthly GIMMS LAI3g GeoTIFF file.
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
% -------------------------------------------------------------------------


%% Step 1: Normalize input values for month and half-month
if isnumeric(month)
    month_str = {'jan', 'feb', 'mar', 'apr', 'may', 'jun', ...
                 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'};
    month = month_str{month};
end

if isnumeric(half_month)
    if half_month == 1, half_month = 'a'; elseif half_month == 2, half_month = 'b'; end
end


%% Step 2: Construct file name: 'AVHRRBUVI04.<year><month><half>.abl'
file_prefix = 'AVHRRBUVI04.';
file_name   = sprintf('%s%d%s%s.abl', file_prefix, year, month, half_month);
file_path   = fullfile(dir_in, file_name);


%% Step 3: Read binary data
fid = fopen(file_path, 'r');
if fid == -1
    error('File %s does not exist.', file_path);
end
data = fread(fid, [2160, 4320], 'uint16=>single', 0, 'ieee-be');
fclose(fid);


%% Step 4: Decode LAI & QC
data(data > 20000) = NaN;

% Get quality flag
qc = data - floor(data/10)*10;
qc(isnan(data)) = NaN;

% Get LAI (scaled value)
lai = floor(data/10)*10.*0.001;

end