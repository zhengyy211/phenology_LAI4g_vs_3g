%% ========================================================================
%  File:    code_3_2_pheno_extract.m
%  Author:  Yaoyao Zheng
%  Email:   zhengyaoyao@stu.pku.edu.cn
%  MATLAB:  R2026a
%
%  Description:
%  Threshold-based extraction of key phenological transition dates. 
%    1) Selection of single-growing-season pixels
%    2) Shifting of VI data for non-summer-green vegetation
%    3) Phenology extraction using threshold methods
% =========================================================================

clear;

%% ------ STEP 1: Load configuration ------
run('code_0_0_config.m');

dir_prep    = cfg.path.prep;
dir_pheno   = cfg.path.pheno;

% List of VI products to process
products    = cfg.vi.products;

% Phenology extraction constraints
GS_min      = cfg.pheno.GS_min;        % Minimum number of growing seasons
GS_max      = cfg.pheno.GS_max;        % Maximum number of growing seasons
buffer_size = cfg.pheno.buffer_size;   % Temporal buffer
trs_list    = cfg.pheno.trs_list;      % Threshold set for phenology extraction
phase_lvls  = cfg.pheno.phase_lvls;    % Phenological phases
phase_lbls  = cfg.pheno.phase_lbls;    % Phenological phases

%% ------ STEP 2: Process each VI product ------
for prod = products
    
    % prod = products(1);
    fprintf('\n[Info] ====== Processing %s ======\n', prod);
    start_time = tic;

    % Time range for processing
    year_start = cfg.vi.(prod).year_start;
    year_end   = cfg.vi.(prod).year_end;
    year_list  = year_start:1:year_end;
    doy_list   = cfg.vi.(prod).doy_list;
    n_years    = numel(year_list);
    n_doys     = numel(doy_list);

    % Savitzky–Golay filter parameters
    if n_doys <= 24
        sg_framelen = 5;
    elseif n_doys <= 46
        sg_framelen = 9;
    else
        sg_framelen = 19;
    end
    
    file_suffix = sprintf('%s_%d_%d', prod, year_start, year_end);

    %% ------ STEP 2.1: Load data ------
    % veg_type : [valid pixel x 1] 1–4 = spring/summer/autumn/winter green
    % data_nGS : [valid pixel x 1] number of growing seasons per year
    % mask_vi  : [all   pixel x 1] logical mask for valid vegetation pixels (copied from VI file)
    file_prior = sprintf('data_pheno_prior_%s.mat', file_suffix);
    s_prior    = load(fullfile(dir_pheno, file_prior), 'data_nGS', 'veg_type');
    % Identify pixels corresponding to a single growing season
    is_single_GS = s_prior.data_nGS > GS_min & s_prior.data_nGS < GS_max;

    file_vi = sprintf('data_%s.mat', file_suffix);
    s_vi    = load(fullfile(dir_prep, file_vi), 'data_vi', 'data_weights', 'mask_vi');
    
    %% ------ STEP 2.2: Select single-growing-season pixels ------ 
    % Update mask to include only pixels with a single growing season
    mask_pheno = s_vi.mask_vi;
    mask_pheno(mask_pheno) = is_single_GS;
    n_pixels = sum(mask_pheno);

    % Get summer-green vegetation mask and corresponding data
    veg_type     = s_prior.veg_type(is_single_GS);
    data_vi      = s_vi.data_vi(is_single_GS, :);
    data_weights = s_vi.data_weights(is_single_GS, :);

    clearvars s_vi

    %% ------ STEP 2.3: Shift VI data according to vegetation type ------
    obs_per_month = n_doys / 12;  % Determine shift (in months) for each vegetation type

    data_vi_shifted      = apply_seasonal_shift(data_vi,      veg_type, obs_per_month);
    data_weights_shifted = apply_seasonal_shift(data_weights, veg_type, obs_per_month);
    
    data_vi_shifted      = num2cell(data_vi_shifted,     2);
    data_weights_shifted = num2cell(data_weights_shifted, 2);

    clearvars data_vi data_weights

    %% ------ STEP 2.4: Extract phenology using threshold method ------ 
    cell_pheno = cell(n_pixels, 1);
    if isempty(gcp('nocreate')), parpool; end
    parfor pixel = 1:n_pixels
        cell_pheno{pixel} = extract_pheno( ...
            data_vi_shifted{pixel}, data_weights_shifted{pixel}, ...
            trs_list, n_years, doy_list, sg_framelen, buffer_size);
    end
    delete(gcp('nocreate')); 

    %% ------ STEP 2.5: Organize extracted phenology data ------ 
    s_pheno = struct;
    for phase = phase_lbls
        s_pheno.(phase) = cell(n_pixels, 1);
    end

    for pixel = 1:n_pixels
        for ii = 1:numel(phase_lvls)
           s_pheno.(phase_lbls(ii)){pixel} = cell_pheno{pixel}.(phase_lvls(ii));
        end
    end

    for phase = phase_lbls
        s_pheno.(phase) = vertcat(s_pheno.(phase){:});
    end

    %% ------ STEP 2.6: Save results ------
    % Save processed data, including:
    % mask_pheno : [all pixel   x 1] logical mask for valid vegetation pixels (True for valid, False for invalid pixels)
    % veg_type   : [valid pixel x 1] 1–4 = spring/summer/autumn/winter green
    % *phases    : [valid pixel x 1] extracted phenological phases (e.g., Greenup, MidGreenup, Peak, MidGreendown, Dormancy, etc.)
    %              These phases represent various key points of the growing season based on VI thresholds
    % n_pixels   : Number of valid pixels with a single growing season
    s_pheno.mask_pheno = mask_pheno;
    s_pheno.veg_type   = s_prior.veg_type(is_single_GS);
    file_pheno         = sprintf('data_pheno_raw_%s.mat', file_suffix);
    save(fullfile(dir_pheno, file_pheno), '-fromstruct', s_pheno, '-v7.3')

    elapsed_time = toc(start_time);
    fprintf('[Info] Processing time for %s: %.2f min\n', prod, elapsed_time / 60);

end