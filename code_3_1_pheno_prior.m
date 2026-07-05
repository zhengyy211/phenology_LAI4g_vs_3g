%% ========================================================================
%  File:    code_3_1_pheno_prior.m
%  Author:  Yaoyao Zheng
%  Email:   zhengyaoyao@stu.pku.edu.cn
%  MATLAB:  R2026a
%
%  Description:
%  Diagnose vegetation growing-season characteristics from VI.
%    1) Classify vegetation growth pattern by dominant green season
%       (spring-green, summer-green, autumn-green, winter-green).
%    2) Estimate number of growing seasons (GS) per year using 
%       smoothed VI peak counts.
% =========================================================================

clear;

%% ------ STEP 1: Load configuration ------
run('code_0_0_config.m');

% Define directories for processed VI data and outputs
dir_prep  = cfg.path.prep;
dir_pheno = cfg.path.pheno;

% List of VI products to process
products  = cfg.vi.products;

% Spatial grid parameters
n_rows    = cfg.space.n_rows;
n_cols    = cfg.space.n_cols;
R         = cfg.space.R;

% Growing season detection parameters
window_size_mth   = cfg.pheno.window_size_mth;    % Moving average window size
min_peak_dist_mth = cfg.pheno.min_peak_dist_mth;  % Minimum distance between peaks for detecting GS


%% ------ STEP 2: Process each VI product ------
for prod = products

    % prod = products(1)
    fprintf('\n[Info] ====== Processing %s ======\n', prod);
    start_time = tic;

    % Time range for processing
    year_start = cfg.vi.(prod).year_start;
    year_end   = cfg.vi.(prod).year_end;
    year_list  = year_start:1:year_end;
    doy_list   = cfg.vi.(prod).doy_list;
    n_years    = numel(year_list);
    n_doys     = numel(doy_list);

    %% ------ STEP 2.1: Load VI data ------
    % Contents of data_<vi_prod>_<year_start>_<year_end>.mat:
    % mask_vi      : [all pixel   x    1] logical mask for valid vegetation pixels
    % data_vi      : [valid pixel x time] filled and corrected VI values
    % data_weights : [valid pixel x time] quality weights for VI values
    % n_pixels     : Number of valid pixels (sum(mask_vi))
    file_suffix = sprintf('%s_%d_%d', prod, year_start, year_end);
    s_vi        = load(fullfile(dir_prep, sprintf('data_%s.mat', file_suffix)));
    n_pixels    = s_vi.n_pixels;

    %% ------ STEP 2.2: Classify vegetation growth type ------
    % Identify the dominant greening season by comparing seasonal mean VI.
    % Categories:
    %   1 = Spring-green
    %   2 = Summer-green
    %   3 = Autumn-green
    %   4 = Winter-green
    % Compute seasonal mean VI
    vi_yearly_mean = reshape(s_vi.data_vi, n_pixels, n_doys, n_years);  % [valid pixel x DOY x year]
    vi_yearly_mean = mean(vi_yearly_mean, 3); % [pixel × DOY]           % [valid pixel x DOY]
    veg_type       = classify_veg_season(vi_yearly_mean, doy_list);

    %% ------ STEP 2.3: Detect number of growing seasons (GS) ------
    % Detect the number of growing seasons (GS) based on peaks
    window_size   = window_size_mth * (n_doys / 12);
    min_peak_dist = min_peak_dist_mth * (n_doys / 12);
    vi_smooth     = movmean(s_vi.data_vi, window_size, 2, 'omitmissing');
    vi_smooth     = num2cell(vi_smooth, 2);
    
    n_peaks = NaN(n_pixels, 1, 'single');
    parfor pixel = 1:n_pixels
        n_peaks(pixel) = detect_num_peaks(vi_smooth{pixel}, min_peak_dist);
    end
    data_nGS = single(n_peaks ./ n_years);  % [valid pixel x 1]

    %% ------ STEP 2.4: Save phenological prior data ------
    % Save processed data, including:
    % veg_type : [valid pixel x 1] 1–4 = spring/summer/autumn/winter green
    % data_nGS : [valid pixel x 1] number of growing seasons per year
    % mask_vi  : [all   pixel x 1] logical mask for valid vegetation pixels (copied from VI file)
    s_out = struct( ...
        'veg_type', veg_type, ...
        'data_nGS', data_nGS, ...
        'mask_vi',  s_vi.mask_vi);

    file_prior = fullfile(dir_pheno, sprintf('data_pheno_prior_%s.mat', file_suffix));
    save(file_prior, '-fromstruct', s_out, '-v7.3')
    
    %% ------ STEP 2.5: Write maps to GeoTIFF ------
    mask_1d = s_vi.mask_vi;

    % Map of vegetation type
    tif_type = fullfile(dir_pheno, sprintf('map_vegetation_type_%s.tif', file_suffix));
    map_type = write_map_to_tiff(veg_type, mask_1d, n_rows, n_cols, R, tif_type);

    % Map of the number of GS per year
    tif_nGS = fullfile(dir_pheno, sprintf('map_nGS_%s.tif', file_suffix));
    map_nGS = write_map_to_tiff(data_nGS, mask_1d, n_rows, n_cols, R, tif_nGS);

    %% ------ STEP 2.6: Plot image for visual check ------
    figure('Name', sprintf('[%s] Vegetation Type', prod), 'Color', 'w');
    h = imagesc(map_type);
    cb = colorbar('southoutside');
    cb.Label.String = 'Dominant Vegetation Type (1–4)';
    cb.FontWeight = 'bold';
    h.AlphaData = ~isnan(map_type);
    title(sprintf('[%s] Dominant Green Season', prod));

    figure('Name', sprintf('[%s] Number of GS', prod), 'Color', 'w');
    h = imagesc(map_nGS);
    cb = colorbar('southoutside');
    cb.Label.String = 'Number of Growing Seasons per Year';
    cb.FontWeight = 'bold';
    h.AlphaData = ~isnan(map_nGS);
    title(sprintf('[%s] Number of GS', prod));

    elapsed_time = toc(start_time);
    fprintf('[Info] Processing time for %s: %.2f min\n', prod, elapsed_time / 60);

end