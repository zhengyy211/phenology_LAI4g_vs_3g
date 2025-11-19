%% ========================================================================
%  File:        code_3_4_pheno_stats.m
%  Author:      Yaoyao Zheng
%  Email:       zhengyaoyao@stu.pku.edu.cn
%  MATLAB:      R2025b
%  Description: This script processes phenology data for various products, 
%               calculates the average phenology values (mean DOY), 
%               computes trends (slopes), and saves the results in spatial maps.
% =========================================================================

clear;

%% ------ STEP 1: Load configuration ------
run('code_0_0_config.m');

% Directories for data
dir_pheno    = cfg.path.pheno;
dir_valid    = cfg.path.valid;

% List of LAI products to process 
% (e.g., ["gimms_lai4g", "gimms_lai3g", "glass_lai", "modis_lai"])
lai_products = cfg.lai.products;
products     = [lai_products, "mcd12q2"];

proc_list    = cfg.proc_list;  % Define product-year combinations

% Spatial / grid info
n_rows       = cfg.space.n_rows;
n_cols       = cfg.space.n_cols;
R            = cfg.space.R;

phase_name_map = containers.Map( ...
    {'SOS_15',  'SOS_50',     'POS',  'EOS_50',       'EOS_15'}, ...
    {'Greenup', 'MidGreenup', 'Peak', 'MidGreendown', 'Dormancy'} ...
);


%% ------ Step 2: Calculate regional average phenology times (DOY) ------
summ_pheno = table();
for prod = products

    % prod = products(1)
    % Time range for processing
    if strcmpi(char(prod), 'modis_lai') | strcmpi(char(prod), 'mcd12q2')
        year_start = cfg.modis.year_start;
        year_end   = cfg.modis.year_end;
    else
        year_start = cfg.lai.year_start;
        year_end   = cfg.lai.year_end;
    end
    year_list = year_start:year_end;
    n_years   = numel(year_list);

    if (ismember(prod, lai_products))
        phases = cfg.pheno.lai.phases;
    else
        phases = cfg.pheno.modis.phases;
    end

    %% Load phenology data
    file_pheno = sprintf('data_pheno_clean_%s_%d_%d.mat', prod, year_start, year_end);
    s_pheno    = load(fullfile(dir_pheno, file_pheno));
    
    for phase = phases

        data_pheno = s_pheno.(phase); % [n_pixels x n_years]

        % ------ Compute Mean DOY (for each pixel over the years) ------
        % Reverse seasonal shifts and calculate the circular mean of DOY
        data_pheno_reverse = apply_seasonal_shift_to_doy(data_pheno, s_pheno.veg_type);
        mean_all = cellfun(@(y) circ_mean_doy(y), num2cell(data_pheno_reverse, 1), 'UniformOutput', false);
        mean_all = round(horzcat(mean_all{:}));

        % Use mapped phase names for output
        if ismember(prod, lai_products)
            name_phase = phase_name_map(phase);
        else
            name_phase = phase;
        end

        temp_summ          = table;
        temp_summ.product  = repmat({prod},       n_years, 1);
        temp_summ.phase    = repmat({name_phase}, n_years, 1);
        temp_summ.year     = year_list(:);
        temp_summ.mean_doy = mean_all(:);

        summ_pheno = [summ_pheno; temp_summ];

    end

end

file_out = fullfile(dir_pheno, 'mean_pheno_all_products.csv');
writetable(summ_pheno, file_out);


%% ------ Step 3: Calculate spatial mean and trend maps ------
for idx = 1:size(proc_list, 1)
    prod          = proc_list{idx, 1};
    yr_start_anal = proc_list{idx, 2};
    yr_end_anal   = proc_list{idx, 3};

    fprintf('\n[Info] ====== Processing %s (%d-%d) ======\n', prod, yr_start_anal, yr_end_anal);
    start_time = tic;

    %% ------ STEP 3.1: Define Data Years for the Analysis ------
    
    
    % Time range for processing
    if strcmpi(char(prod), 'modis_lai') | strcmpi(char(prod), 'mcd12q2')
        yr_start_data = cfg.modis.year_start;
        yr_end_data   = cfg.modis.year_end;
    else
        yr_start_data = cfg.lai.year_start;
        yr_end_data   = cfg.lai.year_end;
    end
    yr_list_data = yr_start_data:yr_end_data;

    if (ismember(prod, lai_products))
        phases = cfg.pheno.lai.phases;
    else
        phases = cfg.pheno.modis.phases;
    end
     
    % Mask for the years of interest
    yr_mask_anal = (yr_list_data >= yr_start_anal) & (yr_list_data <= yr_end_anal);

    file_suffix = sprintf('_%s_%d_%d', prod, yr_start_anal, yr_end_anal);

    %% ------ STEP 3.2: Load phenology data ------
    file_pheno = sprintf('data_pheno_clean_%s_%d_%d.mat', prod, yr_start_data, yr_end_data);
    s_pheno    = load(fullfile(dir_pheno, file_pheno));
    mask_1d    = s_pheno.mask_pheno;


    %% ------ STEP 3.3: Compute and save mean & trend phenology maps ------
    for phase = phases
        
        % Extract the phenology data for the phase (filtered by year range)
        data_pheno = s_pheno.(phase)(:, yr_mask_anal); % [n_pixels x n_years]

        % ------ Mean ------
        mean_all = mean(data_pheno, 2, 'omitmissing');
        mean_all_reverse = apply_seasonal_shift_to_doy(mean_all, s_pheno.veg_type);

        file_mean = fullfile(dir_valid, ['map_pheno_mean_' char(phase) file_suffix '.tif']);
        map_pheno_mean = write_map_to_tiff(mean_all_reverse, mask_1d, n_rows, n_cols, R, file_mean);

        % ------ Trend ------
        % Calculate the linear trend (slope) for each pixel
        [slope_all, pval_all] = calc_linear_trend(data_pheno);

        file_slope = fullfile(dir_valid, ['map_pheno_trend_slope_' char(phase) file_suffix '.tif']);
        file_pval  = fullfile(dir_valid, ['map_pheno_trend_pval_'  char(phase) file_suffix '.tif']);
        write_map_to_tiff(slope_all, mask_1d, n_rows, n_cols, R, file_slope);
        write_map_to_tiff(pval_all,  mask_1d, n_rows, n_cols, R, file_pval);

        % % Plot image for visual check
        % figure('Name', sprintf('[%s] %s', prod, phase), 'Color', 'w');
        % h = imagesc(map_pheno_mean);
        % cb = colorbar('southoutside');
        % cb.Label.String = 'DOY';
        % cb.FontWeight = 'bold';
        % h.AlphaData = ~isnan(map_pheno_mean);
        % title(sprintf('[%s] %s', prod, phase));

    end
    elapsed_time = toc(start_time);
    fprintf('[Info] Processing time: %.2f min\n', elapsed_time/60);

end