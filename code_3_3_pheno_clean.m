%% ========================================================================
%  File:        code_3_3_pheno_clean.m
%  Author:      Yaoyao Zheng
%  Email:       zhengyaoyao@stu.pku.edu.cn
%  MATLAB:      R2025b
%  Description: Postprocesses raw phenology data by detecting outliers.
% =========================================================================

clear;

%% ------ STEP 1: Load configuration ------
run('code_0_0_config.m');

% Define directories for thermal growing season and processed phenology data
dir_prep      = cfg.path.prep;
dir_pheno     = cfg.path.pheno;

% List of LAI products to process 
% (e.g., ["gimms_lai4g", "gimms_lai3g", "glass_lai", "yuan_lai"])
lai_products = cfg.lai.products;
products     = [lai_products, "mcd12q2"];

% Threshold for outlier detection
z_score_thr = cfg.pheno.z_score_thr;

%% ------ STEP 2: Process each product ------
for prod = products

    % prod = products(1)
    fprintf('\n[Info] ====== Processing %s ======\n', prod);
    start_time = tic;

    % Time range for processing
    if strcmpi(char(prod), 'yuan_lai') | strcmpi(char(prod), 'mcd12q2')
        year_start = cfg.modis.year_start;
        year_end   = cfg.modis.year_end;
    else
        year_start = cfg.lai.year_start;
        year_end   = cfg.lai.year_end;
    end

    if (ismember(prod, lai_products))
        phases = cfg.pheno.lai.phases;
    else
        phases = cfg.pheno.modis.phases;
    end

    %% ------ STEP 2.1: Load raw phenology data ------ 
    file_suffix = sprintf('_%s_%d_%d.mat', prod, year_start, year_end);
    s_raw       = load(fullfile(dir_pheno, ['data_pheno_raw' file_suffix]));

    % Initialize the cleaned phenology structure
    s_clean            = struct;
    s_clean.mask_pheno = s_raw.mask_pheno;
    s_clean.veg_type   = s_raw.veg_type;
    n_pixels           = sum(s_clean.mask_pheno);  % Number of valid pixels in the data
    
    %% ------ STEP 2.3: Detect and remove outliers ------
    for phase = phases

        % phase = phases(1)
        data_pheno = s_raw.(phase);
   
        % Detect and Remove Outliers using Z-scores
        % Compute mean and standard deviation per pixel across years
        mean_val = mean(data_pheno, 2, 'omitmissing');
        std_val  =  std(data_pheno, 0, 2, 'omitmissing');

        % Compute Z-scores for outlier detection
        z_scores = abs((data_pheno - mean_val) ./ std_val);
        outliers = z_scores > z_score_thr;

        data_pheno(outliers) = NaN;
        s_clean.(phase) = data_pheno;

    end

    %% ------ STEP 3.4: Save cleaned phenology data ------
    save(fullfile(dir_pheno, ['data_pheno_clean' file_suffix]), '-fromstruct', s_clean, '-v7.3')

    elapsed_time = toc(start_time);
    fprintf('[Info] Processing time: %.2f min\n', elapsed_time / 60);

end