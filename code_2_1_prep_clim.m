%% ========================================================================
%  File:    code_2_1_prep_clim.m
%  Author:  Yaoyao Zheng
%  Email:   zhengyaoyao@stu.pku.edu.cn
%  MATLAB:  R2026a
%  
%  Description:
%  Processes and prepares climate datasets, including CRU TS and ERA5 Land 
%  surface radiation data. Tasks include:
%    1) Loads and preprocesses CRU temperature and precipitation data.
%    2) Performs spatial interpolation and saves the data.
%    3) Processes ERA5 land surface radiation data.
%    4) Generate a thermal growing season mask based on temperature climatology.
% =========================================================================

clear;

%% ------ STEP 1: Load configuration ------
run('code_0_0_config.m');

% Define directories for raw and processed climate data
dir_cru       = fullfile(cfg.path.raw, 'CRU_TS');
dir_era5      = fullfile(cfg.path.raw, 'ERA5_LAND_MONTHLY_AGGR');
dir_prep      = cfg.path.prep;

% Define climatology period (starting one year before analysis period)
yr_st_clim    = cfg.clim.year_start;
yr_ed_clim    = cfg.clim.year_end;
yr_list_clim  = yr_st_clim:1:yr_ed_clim;
n_yrs_clim    = numel(yr_list_clim);
time_suffix   = sprintf('%d_%d', yr_st_clim, yr_ed_clim);

doy_list_clim = cfg.clim.doy_list;
n_doys_clim   = numel(doy_list_clim);

% Spatial grid parameters
n_rows        = cfg.space.n_rows;
n_cols        = cfg.space.n_cols;
lat_limit     = cfg.space.lat_limit;


%% ------ STEP 2: Process CRU TS data ------
% Define the climate variables to process (temperature and precipitation)
vars_clim    = ["temp", "prcp"];
yr_st_cru    = 1901;  % Start year for CRU dataset
yr_ed_cru    = 2024;  % End year for CRU dataset

% Row indices for latitude range
cru_res      = 0.5;  % CRU dataset has 0.5° resolution
cru_lat_vals = ((90 - cru_res/2):-cru_res:(-90 + cru_res/2));
cru_row_st   = find(cru_lat_vals <= lat_limit(2), 1, 'first');
cru_row_ed   = find(cru_lat_vals >= lat_limit(1), 1,  'last');

for clim = vars_clim

    % clim = vars_clim(1);
    fprintf('\n[Info] ====== Processing %s ======\n', clim);
    start_time = tic;

    if clim == "prcp"
        name_clim = "pre";
    elseif clim == "temp"
        name_clim = "tmp";
    end

    %% ------ STEP 2.1: Read CRU NetCDF data ------
    file_cru = sprintf('cru_ts4.09.%d.%d.%s.dat.nc', yr_st_cru, yr_ed_cru, name_clim);
    data_cru = single(ncread(fullfile(dir_cru, file_cru), name_clim));
    
    % Subset the data to match the climatology period
    idx_st   = 12*(yr_st_clim - yr_st_cru) + 1;
    idx_ed   = size(data_cru, 3) - 12*(yr_ed_cru-yr_ed_clim);
    data_cru = rot90(data_cru(:, :, idx_st:idx_ed));
    
    % Restrict to the study domain latitudes (≥ 30°N)
    data_cru = data_cru(cru_row_st:cru_row_ed, :, :);
    
    %% ------ STEP 2.2: Perform Spatial Interpolation ------
    [r, c, t] = size(data_cru, [1, 2, 3]);
    if r == n_rows && c == n_cols
        % do nothing
    elseif r < n_rows && c < n_cols
        cell_cru = cell(t, 1);
        for i_t = 1:t
            temp_cru = imresize(data_cru(:, :, i_t), [n_rows, n_cols], 'bilinear');
            cell_cru{i_t} = reshape(temp_cru, n_rows*n_cols, 1);
        end
        data_cru = horzcat(cell_cru{:});
    elseif r > n_rows && c > n_cols
        data_cru = imresize_block_mean(data_cru, n_rows, n_cols);
        data_cru = reshape(data_cru, n_rows*n_cols, t);
    else
        error('Check data size mismatch.');
    end

    %% ------ STEP 2.3: Save preprocessed CRU data ------
    % Save processed data, including:
    % data : [pixel x time] processed climate data for temperature or precipitation
    %        (interpolated to the target grid and within the climatology period)
    s_out     = struct('data', data_cru);
    file_name = sprintf('data_clim_%s_%s.mat', clim, time_suffix);
    save(fullfile(dir_prep, file_name), '-fromstruct', s_out, '-v7.3')

    clearvars data_cru cell_cru temp_cru

    elapsed_time = toc(start_time);
    fprintf('[Info] Processing time for %s: %.2f min\n', clim, elapsed_time / 60);

end


%% ------ STEP 3: Process ERA5-Land surface radiation ------
fprintf('\n[Info] ====== Processing srad ======\n');
start_time = tic;

% Row indices for latitude range
era_res      = 0.1;  % ERA5-Land dataset has 0.1° resolution
era_lat_vals = ((90 - era_res/2):-era_res:(-90 + era_res/2));
era_row_st   = find(era_lat_vals <= lat_limit(2), 1, 'first');
era_row_ed   = find(era_lat_vals >= lat_limit(1), 1,  'last');

%% ------ STEP 3.1: Read ERA5-Land data ------
% Process ERA5 Land surface radiation data
data_srad = cell(n_yrs_clim, 1);
for idx_yr = 1:n_yrs_clim
    yr = yr_list_clim(idx_yr);
    fprintf('    Year %d...\n', yr);

    % Read the surface radiation data for the current year
    file_era5 = sprintf('ERA5_LAND_MONTHLY_AGGR_srad_%d.tif', yr);
    temp_srad = single(imread(fullfile(dir_era5, file_era5)));

    % Restrict to the study domain latitudes (≥ 30°N)
    data_srad{idx_yr} = temp_srad(era_row_st:era_row_ed, :, :);
end
data_srad = cat(3, data_srad{:});

%% ------ STEP 3.2: Perform Spatial Interpolation ------
[r, c, t] = size(data_srad, [1, 2, 3]);
if r == n_rows && c == n_cols
    % do nothing
elseif r < n_rows && c < n_cols
    cell_srad = cell(t, 1);
    for i_t = 1:t
        temp_srad = imresize(data_srad(:, :, i_t), [n_rows, n_cols], 'bilinear');
        cell_srad{i_t} = reshape(temp_srad, n_rows*n_cols, 1);
    end
    data_srad = horzcat(cell_srad{:});
elseif r > n_rows && c > n_cols
    data_srad = imresize_block_mean(data_srad, n_rows, n_cols);
    data_srad = reshape(data_srad, n_rows*n_cols, t);
else
    error('Check data size mismatch.');
end

% Save processed data, including:
% data : [pixel x time] processed surface radiation data 
%        (interpolated to the target grid)
s_out = struct('data', data_srad);
file_name = sprintf('data_clim_srad_%s.mat', time_suffix);
save(fullfile(dir_prep, file_name), '-fromstruct', s_out, '-v7.3')

clearvars data_srad cell_srad temp_srad

elapsed_time = toc(start_time);
fprintf('[Info] Processing time for srad: %.2f min\n', elapsed_time / 60);

%% ------ STEP 4: Thermal growing season ------
% Load temperature climatology data
s_tmp = load(fullfile(dir_prep, sprintf('data_clim_temp_%s.mat', time_suffix)));

% List of VI products to process
products  = [cfg.vi.products, 'modis_lai'];
for prod = products

    % Time range for processing
    year_start     = cfg.vi.(prod).year_start;
    year_end       = cfg.vi.(prod).year_end;
    year_list      = year_start:1:year_end;
    doy_list       = cfg.vi.(prod).doy_list;
    n_years        = numel(year_list);
    n_doys         = numel(doy_list);

    % Select years that fall within the LAI period
    mask_year      = (yr_list_clim >= year_start) & (yr_list_clim <= year_end);
    mask_mth       = repelem(mask_year, numel(doy_list_clim));
    data_temp      = s_tmp.data(:, mask_mth);
    n_pixels_all   = size(data_temp, 1);

    % Convert valid temperature pixel rows to cell arrays for parallel processing
    is_valid_tmp   = any(~isnan(data_temp), 2);
    n_pixels_valid = sum(is_valid_tmp);
    data_temp      = num2cell(data_temp(is_valid_tmp, :), 2);

    % Construct time axes for interpolation
    offset = (1:n_years)*365 - 365;
    x      = repmat(doy_list_clim, 1, n_years) + repelem(offset, n_doys_clim);
    xq     = repmat(doy_list,      1, n_years) + repelem(offset, n_doys);

    % Perform spline interpolation for temperature data and
    % determine thermal growing season
    thermal_gs_valid = cell(n_pixels_valid, 1);

    if isempty(gcp('nocreate')), parpool; end
    parfor pixel = 1:n_pixels_valid
        temp_interp = interp1(x, data_temp{pixel}, xq, 'spline');
        thermal_gs_valid{pixel} = temp_interp > 0;
    end
    delete(gcp('nocreate')); 

    % Create the final mask for all pixels
    thermal_gs = false(n_pixels_all, numel(xq));
    thermal_gs(is_valid_tmp, :) = vertcat(thermal_gs_valid{:});

    % Save processed data, including:
    % thermal_GS : [all pixel x time] logical mask for temperature > 0°C (thermal growing season)
    s_out     = struct('thermal_gs', thermal_gs);
    file_name = sprintf('data_thermal_gs_%s_%d_%d.mat', prod, year_start, year_end);
    save(fullfile(dir_prep, file_name), '-fromstruct', s_out, '-v7.3')
    
end