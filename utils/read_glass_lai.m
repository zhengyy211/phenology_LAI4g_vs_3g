function lai = read_glass_lai(dir_in, year, doy)
% -------------------------------------------------------------------------
% Function:    read_glass_lai
% Author:      Yaoyao Zheng
% Email:       zhengyaoyao@stu.pku.edu.cn
% Description: This function reads and processes the GLASS LAI data for 
%              a specified year and day-of-year (DOY). It handles missing 
%              or invalid data and scales the LAI values. 
%
% INPUTS:
%   dir_in  - Directory containing the GLASS LAI data files.
%   year    - Year of the data (e.g., 2001).
%   doy     - Day of Year (DOY) (e.g., 211).
%
% OUTPUT:
%   lai     - Processed LAI data (scaled, missing values handled).
% -------------------------------------------------------------------------

%% Step 1: Construct file name and search for matching file
file_prefix = 'GLASS01B02.V40.A';
pattern     = sprintf('%s%d%03d', file_prefix, year, doy);

files = dir(fullfile(dir_in, sprintf('*%s*.hdf', pattern)));
if isempty(files)
    error('No GLASS LAI file found for pattern: %s', pattern);
elseif numel(files) > 1
    warning('Multiple files match pattern "%s", using the first one.', pattern);
end
file_path = fullfile(files(1).folder, files(1).name);


%% Step 2: Read metadata
info  = hdfinfo(file_path);
attrs = struct2table(info.Vgroup.Vgroup(1).SDS.Attributes);

% Extract relevant attributes for scaling and fill value
fill_val     = single(attrs.Value{strcmp(attrs.Name, '_FillValue')});
valid_range  = single(attrs.Value{strcmp(attrs.Name, 'valid_range')});
scale_factor = single(attrs.Value{strcmp(attrs.Name, 'scale_factor')});
add_offset   = single(attrs.Value{strcmp(attrs.Name, 'add_offset')});


%% Step 3: Read LAI data and handle missing/invalid values
lai = single(hdfread(file_path, 'LAI'));

% Replace fill values with NaN
lai(lai == fill_val) = NaN;

% Remove values outside the valid range (set them to NaN)
lai(lai < valid_range(1) | lai > valid_range(2)) = NaN;

% Apply scaling and offset to the LAI data
lai = scale_factor .* lai + add_offset;

end