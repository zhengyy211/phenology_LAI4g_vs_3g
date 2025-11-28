%% ========================================================================
%  File:        code_3_2_pheno_extract.m
%  Author:      Yaoyao Zheng
%  Email:       zhengyaoyao@stu.pku.edu.cn
%  MATLAB:      R2025b
%  Description: Threshold-based extraction of key phenological transition 
%               dates from LAI datasets. 
%               1) Selection of single-growing-season pixels
%               2) Shifting of LAI data for non-summer-green vegetation
%               3) Phenology extraction using threshold methods
% =========================================================================

clear;

%% ------ STEP 1: Load configuration ------
run('code_0_0_config.m');

dir_prep     = cfg.path.prep;
dir_pheno    = cfg.path.pheno;

% List of LAI products to process 
% (e.g., ["gimms_lai4g", "gimms_lai3g", "glass_lai", "yuan_lai"])
lai_products = cfg.lai.products;

% Phenology extraction constraints
GS_min       = cfg.pheno.GS_min;        % Minimum number of growing seasons
GS_max       = cfg.pheno.GS_max;        % Maximum number of growing seasons
sg_framelen  = cfg.pheno.sg_framelen;   % Savitzky–Golay filter window length
buffer_size  = cfg.pheno.buffer_size;   % Temporal buffer
trs_list     = cfg.pheno.lai.trs_list;  % Threshold set for phenology extraction
phases       = cfg.pheno.lai.phases;    % Phenological phases

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

    file_suffix = sprintf('_%s_%d_%d.mat', lai_prod, year_start, year_end);

    %% ------ STEP 2.1: Load data ------
    % veg_type : [valid pixel x 1] 1–4 = spring/summer/autumn/winter green
    % data_nGS : [valid pixel x 1] number of growing seasons per year
    % mask_lai : [all   pixel x 1] logical mask for valid vegetation pixels (copied from LAI file)
    s_prior = load(fullfile(dir_pheno, ['data_pheno_prior' file_suffix]), 'data_nGS', 'veg_type');
    % Identify pixels corresponding to a single growing season
    is_single_GS = s_prior.data_nGS > GS_min & s_prior.data_nGS < GS_max;

    s_lai = load(fullfile(dir_prep, ['data' file_suffix]), 'data_lai',  'data_weights', 'mask_lai');
    
    %% ------ STEP 2.2: Select single-growing-season pixels ------ 
    % Update mask to include only pixels with a single growing season
    mask_pheno = s_lai.mask_lai;
    mask_pheno(mask_pheno) = is_single_GS;
    n_pixels = sum(mask_pheno);

    % Get summer-green vegetation mask and corresponding data
    veg_type     = s_prior.veg_type(is_single_GS);
    data_lai     = s_lai.data_lai(is_single_GS, :);
    data_weights = s_lai.data_weights(is_single_GS, :);

    clearvars s_lai

    %% ------ STEP 2.3: Shift LAI data according to vegetation type ------
    obs_per_month = n_doys / 12;  % Determine shift (in months) for each vegetation type

    data_lai_shifted     = apply_seasonal_shift(data_lai,     veg_type, obs_per_month);
    data_weights_shifted = apply_seasonal_shift(data_weights, veg_type, obs_per_month);
    
    data_lai_shifted     = num2cell(data_lai_shifted,     2);
    data_weights_shifted = num2cell(data_weights_shifted, 2);

    clearvars data_lai data_weights

    %% ------ STEP 2.4: Extract phenology using threshold method ------ 
    cell_pheno = cell(n_pixels, 1);
    parfor pixel = 1:n_pixels
        cell_pheno{pixel} = extract_pheno( ...
            data_lai_shifted{pixel}, data_weights_shifted{pixel}, ...
            trs_list, n_years, doy_list, sg_framelen, buffer_size);
    end

    %% ------ STEP 2.5: Organize extracted phenology data ------ 
    s_pheno = struct;
    for phase = phases
        s_pheno.(phase) = cell(n_pixels, 1);
    end

    for pixel = 1:n_pixels
        for phase = phases
           s_pheno.(phase){pixel} = cell_pheno{pixel}.(phase);
        end
    end

    for phase = phases
        s_pheno.(phase) = vertcat(s_pheno.(phase){:});
    end

    %% ------ STEP 2.6: Save results ------
    % Save processed data, including:
    % mask_pheno : [all pixel   x 1] logical mask for valid vegetation pixels (True for valid, False for invalid pixels)
    % veg_type   : [valid pixel x 1] 1–4 = spring/summer/autumn/winter green
    % *phases    : [valid pixel x 1] extracted phenological phases (e.g., SOS_10, SOS_50, POS, EOS_10, EOS_50, etc.)
    %              These phases represent various key points of the growing season based on LAI thresholds
    % n_pixels   : Number of valid pixels with a single growing season
    s_pheno.mask_pheno = mask_pheno;
    s_pheno.veg_type   = s_prior.veg_type(is_single_GS);
    save(fullfile(dir_pheno, ['data_pheno_raw' file_suffix]), '-fromstruct', s_pheno, '-v7.3')

    elapsed_time = toc(start_time);
    fprintf('[Info] Processing time for %s: %.2f min\n', lai_prod, elapsed_time / 60);

end