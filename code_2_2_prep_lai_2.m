%% ========================================================================
%  File:        code_2_2_prep_lai_2.m
%  Author:      Yaoyao Zheng
%  Email:       zhengyaoyao@stu.pku.edu.cn
%  MATLAB:      R2025b
%  Description: Process LAI products
%               1) Replaces missing values with zero.
%               2) Detects and corrects spikes using a 3-point moving 
%                  average and a threshold based on 25% of the pixel's 
%                  maximum value. Spikes are replaced with the 
%                  minimum of the "good" data.
%               3) Assigns weights based on Quality Control (QC) flags and
%                  fills for each pixel.
%               4) Saves processed pixel-level LAI, weights, and masks for
%                  subsequent analysis.
% =========================================================================

clear;

%% ------ STEP 1: Load configuration ------
run('code_0_0_config.m');

% Define directories for raw LAI data and processed outputs
dir_raw      = cfg.path.raw;
dir_prep     = cfg.path.prep;

% List of LAI products to process 
% (e.g., ["gimms_lai4g", "gimms_lai3g", "glass_lai", "yuan_lai"])
lai_products = cfg.lai.products;

% Spatial / grid info
n_rows       = cfg.space.n_rows;
n_cols       = cfg.space.n_cols;
R            = cfg.space.R;

%% ------ STEP 2: Process each product ------
for lai_prod = lai_products

    fprintf('\n[Info] ====== Processing %s ======\n', lai_prod);
    start_time = tic;

    % Time range for processing
    if strcmpi(char(lai_prod), 'yuan_lai')
        year_start = cfg.modis.year_start;
        year_end   = cfg.modis.year_end;
    else
        year_start = cfg.lai.year_start;
        year_end   = cfg.lai.year_end;
    end
    year_list = year_start:1:year_end;
    n_years   = numel(year_list);

    % Determine the list of DOYs based on the product
    if strcmpi(char(lai_prod), 'gimms_lai3g') | strcmpi(char(lai_prod), 'gimms_lai4g')
        doy_list = cfg.lai.gimms.doy_list;
    elseif strcmpi(char(lai_prod), 'glass_lai') | strcmpi(char(lai_prod), 'yuan_lai')
        doy_list = cfg.lai.glass_yuan.doy_list;
    end
    n_doys = numel(doy_list);

    %% ------ STEP 2.1: Load Raw LAI Data ------
    % Load raw LAI data along with its QC and mask information
    %   data_lai  : [valid_pixel x time] (rows correspond to true entries of mask_lai)
    %   data_qc   : [valid_pixel x time] quality control flags
    %   mask_lai  : [all_pixel   x    1] logical mask for valid LAI pixels
    file_raw = sprintf('data_raw_%s_%d_%d.mat', lai_prod, year_start, year_end);
    load(fullfile(dir_raw, file_raw), 'data_lai', 'data_qc', 'mask_lai');

    %% ------ STEP 2.2: Load thermal growing season ------
    % The file provides:
    % thermal_GS : [all pixel x time] logical mask for temperature > 0°C (thermal growing season)
    prod_type = strsplit(lai_prod, '_');
    prod_type = prod_type{1};
    file_tmp = sprintf('data_thermal_GS_%s_%d_%d.mat', prod_type, year_start, year_end);

    load(fullfile(dir_prep, file_tmp), 'thermal_GS');
    thermal_GS = thermal_GS(mask_lai, :);

    %% ------ STEP 2.3: Fill missing values ------
    mask_fill_missing = isnan(data_lai);
    data_lai(mask_fill_missing) = 0;

    %% ------ STEP 2.4: Spike detection and correction ------
    % A spike is defined as when the absolute difference between the observed LAI
    % and the 3-point moving average exceeds 25% of the pixel's maximum value.
    good_qc = 0;
    [data_lai, mask_spike_all] = remove_spike(data_lai, data_qc, good_qc, thermal_GS);

    %% ------ STEP 2.5: Assign weights according to QC ------
    data_weights = NaN(size(data_qc), 'single');
    data_weights(data_qc == 0)     =   1;  % good quality
    data_weights(data_qc == 1)     = 0.8;  % spline interpolation
    data_weights(data_qc == 2)     = 0.5;  % possible snow/cloud cover
    data_weights(data_qc == 3)     = 0.2;  % missing

    % Further downweight filled or spike-corrected pixels
    data_weights(mask_fill_missing) = 0.2;
    data_weights(mask_spike_all)    = 0.2;
    data_weights(isnan(data_qc))    = 0.2;

    clearvars mask_fill_missing mask_spike_all

    %% ------ STEP 2.7: Save processed data ------
    % Save processed data, including:
    % mask_lai     : [all pixel   x    1] logical mask for valid vegetation pixels
    % data_lai     : [valid pixel x time] filled and corrected LAI values
    % data_weights : [valid pixel x time] corresponding weights for LAI values
    % n_pixels     : Number of valid pixels
    n_pixels = sum(mask_lai);

    file_lai = sprintf('data_%s_%d_%d.mat', lai_prod, year_start, year_end);
    save(fullfile(dir_prep, file_lai), 'data_lai', 'data_weights', 'mask_lai', 'n_pixels', '-v7.3')

    %% ------ STEP 2.8: Plot regional mean LAI for visual check ------
    mean_lai = mean(data_lai, 1, 'omitmissing');

    figure('Name', sprintf('[%s] Mean LAI', lai_prod), 'Color', 'w');
    plot(mean_lai, '-b', 'LineWidth', 1.5);
    ylabel('LAI (m^2/m^2)');
    title(sprintf('[%s] Regional mean LAI', lai_prod));
    grid on;

    elapsed_time = toc(start_time);
    fprintf('[Info] Processing time for %s: %.2f min\n', lai_prod, elapsed_time / 60);

    clearvars data_lai data_weights mask_lai

end