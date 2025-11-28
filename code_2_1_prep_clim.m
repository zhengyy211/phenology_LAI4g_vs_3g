%% ========================================================================
%  File:        code_2_1_prep_clim.m
%  Author:      Yaoyao Zheng
%  Email:       zhengyaoyao@stu.pku.edu.cn
%  MATLAB:      R2025b
%  Description: Processes and prepares climate datasets, including CRU TS 
%               and ERA5 Land surface radiation data. Tasks include:
%               1) Loads and preprocesses CRU temperature and precipitation data.
%               2) Performs spatial interpolation and saves the data.
%               3) Processes ERA5 land surface radiation data.
%               4) Generate a thermal growing season mask based on temperature 
%                  climatology.
%% ========================================================================

clear;

%% ------ STEP 1: Load configuration ------
run('code_0_0_config.m');

% Define directories for raw and processed climate data
dir_cru  = fullfile(cfg.path.raw, 'cru_ts');
dir_era5 = fullfile(cfg.path.raw, 'era5_land');
dir_prep = cfg.path.prep;

% Define climatology period (starting one year before analysis period)
yr_start_clim  = cfg.clim.year_start;
yr_end_clim    = cfg.clim.year_end;
yr_list_clim   = yr_start_clim:1:yr_end_clim;
n_yrs_clim     = numel(yr_list_clim);
time_suffix    = sprintf('%d_%d', yr_start_clim, yr_end_clim);

doy_list_clim  = cfg.clim.doy_list;
n_doys_clim    = numel(doy_list_clim);

% Spatial grid parameters
n_rows    = cfg.space.n_rows;
n_cols    = cfg.space.n_cols;
lat_limit = cfg.space.lat_limit;


%% ------ STEP 2: Process CRU TS data ------

% Define the climate variables to process (temperature and precipitation)
vars_clim    = ["tmp", "prcp"];
cru_yr_start = 1901;  % Start year for CRU dataset
cru_yr_end   = 2024;  % End year for CRU dataset

% Row indices for latitude range
cru_spatial_res = 0.5;  % CRU dataset has 0.5° resolution
cru_lat_vals    = ((90 - cru_spatial_res/2):-cru_spatial_res:(-90 + cru_spatial_res/2));
cru_row_st      = find(cru_lat_vals <= lat_limit(2), 1, 'first');
cru_row_ed      = find(cru_lat_vals >= lat_limit(1), 1,  'last');

for clim = vars_clim

    fprintf('\n[Info] ====== Processing %s ======\n', clim);

    if clim == "prcp"
        name_clim = "pre";
    elseif clim == "tmp"
        name_clim = clim;
    end

    %% ------ STEP 2.1: Read CRU NetCDF data ------
    file_cru = sprintf('cru_ts4.09.%d.%d.%s.dat.nc', cru_yr_start, cru_yr_end, name_clim);
    data_cru = single(ncread(fullfile(dir_cru, file_cru), name_clim));
    
    % Subset the data to match the climatology period
    idx_st   = 12*(yr_start_clim - cru_yr_start) + 1;
    idx_ed   = size(data_cru, 3) - 12*(cru_yr_end-yr_end_clim);
    data_cru = rot90(data_cru(:, :, idx_st:idx_ed));
    
    % Restrict to the study domain latitudes (≥ 30°N)
    data_cru = data_cru(cru_row_st:cru_row_ed, :, :);
    
    %% ------ STEP 2.2: Perform Spatial Interpolation ------
    % Perform bilinear interpolation to match the target grid size
    n_month = size(data_cru, 3);
    cell_clim = cell(n_month, 1);
    for idx_mth = 1:n_month
        temp_clim = imresize(data_cru(:, :, idx_mth), [n_rows, n_cols], 'bilinear');
        cell_clim{idx_mth} = reshape(temp_clim, [], 1);
    end
    data_clim = horzcat(cell_clim{:});

    %% ------ STEP 2.3: Save preprocessed CRU data ------
    % Save processed data, including:
    % data : [pixel x time] processed climate data for temperature or precipitation
    %        (interpolated to the target grid and within the climatology period)
    s_out = struct('data', data_clim);
    file_name = ['data_clim_' char(clim) '_' time_suffix '.mat'];
    save(fullfile(dir_prep, file_name), '-fromstruct', s_out, '-v7.3')

    clearvars data_cru cell_clim data_clim temp_clim

end

%% ------ STEP 3: Process ERA5 Land surface radiation ------
fprintf('\n[Info] ====== Processing srad ======\n');

% Process ERA5 Land surface radiation data
cell_srad = cell(n_yrs_clim, 1);
for idx_yr = 1:n_yrs_clim
    current_year = yr_list_clim(idx_yr);
    fprintf('[Info] Year %d...\n', current_year);

    % Read the surface radiation data for the current year
    file_era5 = sprintf('era5_land_srad_%04d.tif', current_year);
    data_srad = single(imread(fullfile(dir_era5, file_era5)));
    cell_srad{idx_yr} = reshape(data_srad, [], size(data_srad, 3));
end
data_clim = horzcat(cell_srad{:});

% Save processed data, including:
% data : [pixel x time] processed surface radiation data 
%        (interpolated to the target grid)
s_out = struct('data', data_clim);
file_name = ['data_clim_srad_' time_suffix '.mat'];
save(fullfile(dir_prep, file_name), '-fromstruct', s_out, '-v7.3')

clearvars cell_srad data_srad data_clim


%% ------ STEP 4: Thermal growing season ------
% Load temperature climatology data
s_tmp = load(fullfile(dir_prep, ['data_clim_tmp_' time_suffix '.mat']));

for prod = ["gimms", "glass", "yuan"]

    % Time range for processing
    if strcmpi(char(prod), 'yuan')
        yr_start_lai = cfg.modis.year_start;
        yr_end_lai   = cfg.modis.year_end;
    else
        yr_start_lai = cfg.lai.year_start;
        yr_end_lai   = cfg.lai.year_end;
    end
    n_yrs_lai = yr_end_lai - yr_start_lai + 1;

    % Determine the list of DOYs based on the product
    if strcmpi(char(prod), 'gimms')
        doy_list_lai = cfg.lai.gimms.doy_list;
    elseif strcmpi(char(prod), 'glass') | strcmpi(char(prod), 'yuan')
        doy_list_lai = cfg.lai.glass_yuan.doy_list;
    end
    n_doys_lai = numel(doy_list_lai);

    % Select years that fall within the LAI period
    mask_year     = (yr_list_clim >= yr_start_lai) & (yr_list_clim <= yr_end_lai);
    mask_mth      = repelem(mask_year, numel(doy_list_clim));
    data_tmp      = s_tmp.data(:, mask_mth);
    n_pixels_all  = size(data_tmp, 1);

    % Convert valid temperature pixel rows to cell arrays for parallel processing
    is_valid_tmp   = any(~isnan(data_tmp), 2);
    n_pixels_valid = sum(is_valid_tmp);
    data_tmp       = num2cell(data_tmp(is_valid_tmp, :), 2);

    % Construct time axes for interpolation
    offset = (1:n_yrs_lai)*365 - 365;
    x  = repmat(doy_list_clim, 1, n_yrs_lai) + repelem(offset, n_doys_clim);
    xq = repmat(doy_list_lai,  1, n_yrs_lai) + repelem(offset, n_doys_lai);

    % Perform spline interpolation for temperature data and
    % determine thermal growing season
    thermal_GS_valid = cell(n_pixels_valid, 1);
    parfor pixel = 1:n_pixels_valid
        temp_interp = interp1(x, data_tmp{pixel}, xq, 'spline');
        thermal_GS_valid{pixel} = temp_interp > 0;
    end

    % Create the final mask for all pixels
    thermal_GS = false(n_pixels_all, numel(xq));
    thermal_GS(is_valid_tmp, :) = vertcat(thermal_GS_valid{:});

    % Save processed data, including:
    % thermal_GS : [all pixel x time] logical mask for temperature > 0°C (thermal growing season)
    s_out = struct('thermal_GS', thermal_GS);
    file_name = sprintf('data_thermal_GS_%s_%d_%d.mat', prod, yr_start_lai, yr_end_lai);
    save(fullfile(dir_prep, file_name), '-fromstruct', s_out, '-v7.3')
end