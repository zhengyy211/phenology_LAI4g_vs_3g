function [lai, qc] = read_gimms_lai3g(dir_in, year, month, half_month)
% -------------------------------------------------------------------------
% Function:    read_gimms_lai3g
% Author:      Yaoyao Zheng
% Email:       zhengyaoyao@stu.pku.edu.cn
% Description: This function reads and processes the GIMMS LAI3g data for 
%              a specified year, month, and half-month period. It handles
%              missing or invalid data and scales the LAI values, and extracts the QC flags.
%
% INPUTS:
%   dir_in      - Directory containing the GIMMS LAI3g data files.
%   year        - Year of the data (e.g., 2001).
%   month       - Month of the data (1-12).
%   half_month  - Half-month period (1 or 2).
%
% OUTPUT:
%   lai         - Processed LAI data (scaled, missing values handled).
%   qc          - Quality control flags for the LAI data.
% -------------------------------------------------------------------------

%% Step 1: Normalize input values for month and half-month
if isnumeric(month)
    month_str = {'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'};
    month = month_str{month};
end

if isnumeric(half_month)
    if half_month == 1
        half_month = 'a';
    elseif half_month == 2
        half_month = 'b';
    end
end


%% Step 2: Construct file name: 'AVHRRBUVI04.<year><month><half>.abl'
file_prefix = 'AVHRRBUVI04.';
file_name   = sprintf('%s%d%s%s.abl', file_prefix, year, month, half_month);
file_path   = fullfile(dir_in, file_name);


%% Step 3: Read binary data
fid = fopen(file_path,'r');
if fid == -1
    error('File %s does not exist.', file_path);
end
data = fread(fid, [2160, 4320], 'uint16=>single', 0, 'ieee-be');
fclose(fid);


%% Step 4: Decode LAI & QC
data(data > 20000) = nan;

% Get quality flag
qc = data - floor(data/10)*10;
qc(isnan(data)) = NaN;

% Get LAI (scaled value)
lai = floor(data/10)*10.*0.001;

end