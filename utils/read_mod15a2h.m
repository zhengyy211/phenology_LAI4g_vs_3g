function lai = read_mod15a2h(dir_in, year)
% -------------------------------------------------------------------------
% Function:    read_mod15a2h
% Author:      Yaoyao Zheng
% Email:       zhengyaoyao@stu.pku.edu.cn
% Description: This function reads and processes the MOD15A2H data for
%              a specified year.
%
% INPUTS:
%   dir_in  - Directory containing the MOD15A2H data files.
%   year    - Year of the data (e.g., 2001).
%
% OUTPUT:
%   lai     - Processed LAI data (missing values handled).
% -------------------------------------------------------------------------

%% Step 1: Construct file name and search for matching file
file_pattern = sprintf('MOD15A2H_%d', year);
file_path_W  = fullfile(dir_in, [file_pattern, '_W.tif']);  % Western tile
file_path_E  = fullfile(dir_in, [file_pattern, '_E.tif']);  % Eastern tile

if ~isfile(file_path_W)
    error('File not found: %s', file_path_W);
end
if ~isfile(file_path_E)
    error('File not found: %s', file_path_E);
end


%% Step 2: Read and merge LAI tiles
lai_W = single(imread(file_path_W));
lai_E = single(imread(file_path_E));

% Merge west and east tiles horizontally
lai = cat(2, lai_W, lai_E);   % size: [H, W_total, num_scenes]


%% Step 3: Check scene count (temporal dimension)
n_scenes = size(lai, 3);
if year == 2001 && n_scenes == 44
    warning('Year 2001: Only 44 scenes found. Performing linear interpolation to 46 scenes.');

    doy_full    = 1:8:361;
    doy_missing = [169, 177];
    doy_exist   = setdiff(doy_full, doy_missing);

    lai_full = NaN(size(lai, 1), size(lai, 2), numel(doy_full), 'like', lai);
    lai_full(:,:, ismember(doy_full, doy_exist)) = lai;

    idx161 = doy_full == 161;
    idx169 = doy_full == 169;
    idx177 = doy_full == 177;
    idx185 = doy_full == 185;

    diff185_161 = lai_full(:, :, idx185) - lai_full(:, :, idx161);

    lai_full(:, :, idx169) = lai_full(:, :, idx161) + diff185_161 * (1/3);
    lai_full(:, :, idx177) = lai_full(:, :, idx161) + diff185_161 * (2/3);

    lai = lai_full;
elseif size(lai, 3) ~= 46
    error('Year %d: Expected 46 scenes, but found %d scenes.', year, n_scenes);
end
end