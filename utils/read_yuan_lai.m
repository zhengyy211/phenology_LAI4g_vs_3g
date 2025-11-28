function lai = read_yuan_lai(dir_in, year)
% -------------------------------------------------------------------------
% Function:    read_yuan_lai
% Author:      Yaoyao Zheng
% Email:       zhengyaoyao@stu.pku.edu.cn
% Description: This function reads and processes the Yuan LAI data for 
%              a specified year. It handles missing. 
%
% INPUTS:
%   dir_in  - Directory containing the Yuan LAI data files.
%   year    - Year of the data (e.g., 2001).
%
% OUTPUT:
%   lai     - Processed LAI data (missing values handled).
% -------------------------------------------------------------------------

%% Step 1: Construct file name and search for matching file
file_prefix = 'lai_8-day_0.1_';
file_name   = sprintf('%s%d.nc', file_prefix, year);
file_path   = fullfile(dir_in, file_name);

if ~isfile(file_path)
    error('File %s does not exist.', file_path);
end


%% Step 2: Read metadata
info     = ncinfo(file_path);
attrs    = struct2table(info.Variables);
fill_val = single(attrs.FillValue(strcmp(attrs.Name, 'lai')));


%% Step 3: Read LAI data and handle missing/invalid values
lai = single(flipud(rot90(ncread(file_path, 'lai'))));


% Replace fill values with NaN
lai(lai == fill_val) = NaN;

end