%% ========================================================================
%  File:    code_3_3_pheno_clean.m
%  Author:  Yaoyao Zheng
%  Email:   zhengyaoyao@stu.pku.edu.cn
%  MATLAB:  R2026a
%  
%  Description:
%  Postprocesses raw phenology data by detecting outliers.
% =========================================================================

clear;

%% ------ STEP 1: Load configuration ------
run('code_0_0_config.m');

% Define directories for thermal growing season and processed phenology data
dir_prep      = cfg.path.prep;
dir_pheno     = cfg.path.pheno;

% List of VI products to process
products      = cfg.vi.products;
products      = [products, "mcd12q2"];

year_mid      = cfg.time.year_mid;

phases        = cfg.pheno.phase_lbls;
n_phases      = numel(phases);

% Threshold for outlier detection
z_score_thr   = cfg.pheno.z_score_thr;

%% ------ STEP 2: Process each product ------
for prod = products

    % prod = products(1);
    fprintf('\n[Info] ====== Processing %s ======\n', prod);
    start_time = tic;

    % Time range for processing
    year_start = cfg.vi.(prod).year_start;
    year_end   = cfg.vi.(prod).year_end;
    year_list  = year_start:1:year_end;

    %% ------ STEP 2.1: Load raw phenology data ------ 
    file_suffix = sprintf('%s_%d_%d', prod, year_start, year_end);
    file_raw    = sprintf('data_pheno_raw_%s.mat', file_suffix);
    s_raw       = load(fullfile(dir_pheno, file_raw));

    %% ------ STEP 2.3: Detect and remove outliers ------
    s_clean = struct;
    cell_valid = cell(n_phases, 1);
    for idx_phase = 1:n_phases

        phase = phases(idx_phase);
        data_pheno = s_raw.(phases(idx_phase));
   
        % Detect and Remove Outliers using Z-scores
        % Compute mean and standard deviation per pixel across years
        mean_val = mean(data_pheno,    2, 'omitmissing');
        std_val  =  std(data_pheno, 0, 2, 'omitmissing');

        % Compute Z-scores for outlier detection
        z_scores = abs((data_pheno - mean_val) ./ std_val);
        outliers = z_scores > z_score_thr;

        data_pheno(outliers) = NaN;
        s_clean.(phase) = data_pheno;

        if year_start == cfg.time.year_start
            cell_valid{idx_phase} = sum(~isnan(data_pheno(:, year_list <= year_mid)), 2) >= 10 & ...
                                    sum(~isnan(data_pheno(:, year_list >  year_mid)), 2) >= 10; 
        else
            cell_valid{idx_phase} = sum(~isnan(data_pheno(:, year_list >  year_mid)), 2) >= 10; 
        end

    end

    mask_valid = horzcat(cell_valid{:});
    mask_valid = all(mask_valid, 2);
    for phase = phases
        s_clean.(phase) = s_clean.(phase)(mask_valid, :);
    end

    mask_pheno = s_raw.mask_pheno;
    mask_pheno(mask_pheno) = mask_valid;
    veg_type = s_raw.veg_type(mask_valid, :);

    s_clean.veg_type   = veg_type;
    s_clean.mask_pheno = mask_pheno;

    %% ------ STEP 3.4: Save cleaned phenology data ------
    file_clean = sprintf('data_pheno_clean_%s.mat', file_suffix);
    save(fullfile(dir_pheno, file_clean), '-fromstruct', s_clean, '-v7.3')

    elapsed_time = toc(start_time);
    fprintf('[Info] Processing time: %.2f min\n', elapsed_time / 60);

end