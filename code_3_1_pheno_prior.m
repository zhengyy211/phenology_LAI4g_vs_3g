%% ========================================================================
%  File:        code_3_1_pheno_prior.m
%  Author:      Yaoyao Zheng
%  Email:       zhengyaoyao@stu.pku.edu.cn
%  MATLAB:      R2025b
%  Description: Diagnose vegetation growing-season characteristics from LAI.
%               1) Classify vegetation growth pattern by dominant green season
%                  (spring-green, summer-green, autumn-green, winter-green).
%               2) Estimate number of growing seasons (GS) per year using
%                  smoothed LAI peak counts.
% =========================================================================

clear;

%% ------ STEP 1: Load configuration ------
run('code_0_0_config.m');

% Define directories for processed LAI data and outputs
dir_prep      = cfg.path.prep;
dir_pheno     = cfg.path.pheno;

% List of LAI products to process 
% (e.g., ["gimms_lai4g", "gimms_lai3g", "glass_lai", "yuan_lai"])
lai_products  = cfg.lai.products;

% Spatial grid parameters
n_rows        = cfg.space.n_rows;
n_cols        = cfg.space.n_cols;
R             = cfg.space.R;

% Growing season detection parameters
window_size_mth   = cfg.pheno.window_size_mth;    % Moving average window size
min_peak_dist_mth = cfg.pheno.min_peak_dist_mth;  % Minimum distance between peaks for detecting GS

%% ------ STEP 2: Process each LAI product ------
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

    %% ------ STEP 2.1: Load LAI data ------
    % Contents of data_<lai_prod>_<year_start>_<year_end>.mat:
    % mask_lai     : [all pixel   x    1] logical mask for valid vegetation pixels
    % data_lai     : [valid pixel x time] filled and corrected LAI values
    % data_weights : [valid pixel x time] quality weights for LAI values
    % n_pixels     : Number of valid pixels (sum(mask_lai))
    file_suffix = sprintf('_%s_%d_%d', lai_prod, year_start, year_end);
    s_lai       = load(fullfile(dir_prep, ['data' file_suffix '.mat']));
    n_pixels    = s_lai.n_pixels;


    %% ------ STEP 2.2: Classify vegetation growth type ------
    % Identify the dominant greening season by comparing seasonal mean LAI.
    % Categories:
    %   1 = Spring-green
    %   2 = Summer-green
    %   3 = Autumn-green
    %   4 = Winter-green
    % Compute seasonal mean LAI
    lai_yearly_mean = reshape(s_lai.data_lai, n_pixels, n_doys, n_years);  % [valid pixel x DOY x year]
    lai_yearly_mean = mean(lai_yearly_mean, 3); % [pixel × DOY]            % [valid pixel x DOY]
    veg_type        = classify_veg_season(lai_yearly_mean, doy_list);

    %% ------ STEP 2.3: Detect number of growing seasons (GS) ------
    % Detect the number of growing seasons (GS) based on peaks
    window_size   = window_size_mth * (n_doys / 12);
    min_peak_dist = min_peak_dist_mth * (n_doys / 12);
    lai_smooth    = movmean(s_lai.data_lai, window_size, 2, 'omitmissing');
    lai_smooth    = num2cell(lai_smooth, 2);
    
    n_peaks = NaN(n_pixels, 1, 'single');
    parfor pixel = 1:n_pixels
        n_peaks(pixel) = detect_num_peaks(lai_smooth{pixel}, min_peak_dist);
    end
    data_nGS = single(n_peaks ./ n_years);  % [valid pixel x 1]

    %% ------ STEP 2.4: Save phenological prior data ------
    % Save processed data, including:
    % veg_type : [valid pixel x 1] 1–4 = spring/summer/autumn/winter green
    % data_nGS : [valid pixel x 1] number of growing seasons per year
    % mask_lai : [all   pixel x 1] logical mask for valid vegetation pixels (copied from LAI file)
    s_out = struct( ...
        'veg_type', veg_type, ...
        'data_nGS', data_nGS, ...
        'mask_lai', s_lai.mask_lai);

    file_prior = fullfile(dir_pheno, ['data_pheno_prior' file_suffix '.mat']);
    save(file_prior, '-fromstruct', s_out, '-v7.3')
    
    %% ------ STEP 2.5: Write maps to GeoTIFF ------
    mask_1d = s_lai.mask_lai;

    % Map of vegetation type
    tif_type = fullfile(dir_pheno, ['map_vegetation_type' file_suffix '.tif']);
    map_type = write_map_to_tiff(veg_type, mask_1d, n_rows, n_cols, R, tif_type);

    % Map of the number of GS per year
    tif_nGS = fullfile(dir_pheno, ['map_nGS' file_suffix '.tif']);
    map_nGS = write_map_to_tiff(data_nGS, mask_1d, n_rows, n_cols, R, tif_nGS);

    %% ------ STEP 2.5: Plot image for visual check ------
    figure('Name', sprintf('[%s] Vegetation Type', lai_prod), 'Color', 'w');
    h = imagesc(map_type);
    cb = colorbar('southoutside');
    cb.Label.String = 'Dominant Vegetation Type (1–4)';
    cb.FontWeight = 'bold';
    h.AlphaData = ~isnan(map_type);
    title(sprintf('[%s] Dominant Green Season', lai_prod));

    figure('Name', sprintf('[%s] Number of GS', lai_prod), 'Color', 'w');
    h = imagesc(map_nGS);
    cb = colorbar('southoutside');
    cb.Label.String = 'Number of Growing Seasons per Year';
    cb.FontWeight = 'bold';
    h.AlphaData = ~isnan(map_nGS);
    title(sprintf('[%s] Number of GS', lai_prod));

    elapsed_time = toc(start_time);
    fprintf('[Info] Processing time for %s: %.2f min\n', lai_prod, elapsed_time / 60);

end