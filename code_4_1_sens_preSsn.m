%% ========================================================================
%  File:        code_4_1_sens_preSsn.m
%  Author:      Yaoyao Zheng
%  Email:       zhengyaoyao@stu.pku.edu.cn
%  MATLAB:      R2025b
%  Description: Calculates pre-season climate variables for phenology sensitivity
%               analysis. For each LAI product and phenological phase, the script:
%               1) Extracts cleaned phenology data and corresponding climate data,
%               2) Computes climate time series at multiple pre-season lags,
%               3) Identifies the optimal lag maximizing correlation with phenology,
%               4) Stores both the optimal lag and aligned pre-season climate data.
%               Outputs are saved for subsequent climate sensitivity calculations.
% =========================================================================

clear;

%% ------ STEP 1: Load configuration ------
run('code_0_0_config.m');

% Define directories for pre-processed data, phenology data, and sensitivity results
dir_prep       = cfg.path.prep;
dir_pheno      = cfg.path.pheno;
dir_sens       = cfg.path.sens;

% Define variables for processing
lai_products   = cfg.lai.products;
products       = [lai_products, "mcd12q2"];

vars_clim      = cfg.sens.vars_clim ;  % Climate variables to process
year_st_clim   = cfg.clim.year_start;  % climate data starts 1 year earlier for lag
year_ed_clim   = cfg.clim.year_end;
year_list_clim = year_st_clim:year_ed_clim;

max_lag_mths   = cfg.sens.max_lag_mths;  % Maximum number of lag periods to consider


%% ------ STEP 2: Loop through each product (LAI or MODIS) ------
for prod = products

    fprintf('\n[Info] ====== Processing %s ======\n', prod);
    start_time = tic;

    %% ------ STEP 2.1: Determine data years ------

    % Time range for processing
    if strcmpi(char(prod), 'yuan_lai') | strcmpi(char(prod), 'mcd12q2')
        year_st_pheno = cfg.modis.year_start;
        year_ed_pheno = cfg.modis.year_end;
    else
        year_st_pheno = cfg.lai.year_start;
        year_ed_pheno = cfg.lai.year_end;
    end
    year_list_pheno = year_st_pheno:year_ed_pheno;
    n_years         = numel(year_list_pheno);

    if (ismember(prod, lai_products))
        phases = cfg.pheno.lai.phases;
    else
        phases = cfg.pheno.modis.phases;
    end

    % Mask for climate years (consider the previous year for lag)
    year_mask_clim  = (year_list_clim >= (year_st_pheno - 1)) & (year_list_clim <= year_ed_pheno);
    month_mask_clim = repelem(year_mask_clim, 12);

    %% ------ STEP 2.2: Load phenology data ------
    file_pheno = sprintf('data_pheno_clean_%s_%d_%d.mat', prod, year_st_pheno, year_ed_pheno);
    s_pheno = load(fullfile(dir_pheno, file_pheno));
    n_pixel = sum(s_pheno.mask_pheno);

    for lag = max_lag_mths

        fprintf('\n[Info] Lag %d\n', lag);

        s_preSsn = struct;
        s_preSsn.mask_pheno = s_pheno.mask_pheno;

        for clim = vars_clim

            if clim == "prcp"
                method = 'sum';
            else
                method = 'mean';
            end

            %% ------ STEP 2.3:  Process climate data for each climate variable ------
            file_clim     = sprintf('data_clim_%s_%d_%d.mat', clim, year_st_clim, year_ed_clim);
            s_clim        = load(fullfile(dir_prep, file_clim));
            data_clim_raw = s_clim.data(s_pheno.mask_pheno, :);
            clearvars s_clim

            % Initialize the empty matrix for shifted data
            obs_per_month     = 1;  % 1 observation per month
            data_clim_shifted = apply_seasonal_shift(data_clim_raw, s_pheno.veg_type, obs_per_month);
            data_clim_shifted = data_clim_shifted(:, month_mask_clim);

            % Identify pixels with valid climate data
            valid_clim = any(~isnan(data_clim_shifted), 2);
            data_clim_shifted = num2cell(data_clim_shifted, 2);

            for phase = phases

                fprintf('[Info] %s - %s...\n', phase, clim);

                % Get phenology data for the current phase and year range
                data_pheno  = s_pheno.(phase);
                valid_pheno = sum(~isnan(data_pheno), 2) > 3;
                data_pheno  = num2cell(data_pheno, 2);

                valid_pixels   = valid_pheno & valid_clim;
                n_valid_pixels = sum(valid_pixels);
                data_pheno     = data_pheno(valid_pixels);
                data_clim      = data_clim_shifted(valid_pixels);

                %% ------ STEP 2.4: Parallel computation of pre-season climate ----
                data_best_lag = NaN(n_valid_pixels, 1, 'single');
                cell_clim_lag = cell(n_valid_pixels, 1);

                parfor i_pixel = 1:n_valid_pixels

                    % Calculate the best pre-season lag using climate-phenology correlation
                    [data_best_lag(i_pixel), cell_clim_lag{i_pixel}] = ...
                        calc_preSSn_clim(data_pheno{i_pixel}, data_clim{i_pixel}, n_years, lag, method);

                end

                %% ------ STEP 2.5: Organize results ------
                data_best_lag_all = NaN(n_pixel, 1, 'single');
                data_best_lag_all(valid_pixels) = data_best_lag;

                data_clim_lag_all = NaN(n_pixel, n_years,'single');
                data_clim_lag_all(valid_pixels, :) = vertcat(cell_clim_lag{:});

                s_preSsn.(phase).(sprintf('best_lag_%s', clim)) = data_best_lag_all;
                s_preSsn.(phase).(sprintf('preSsn_%s',   clim)) = data_clim_lag_all;
                s_preSsn.(phase).('pheno')                      = s_pheno.(phase);

            end

        end

        file_preSsn = sprintf('data_preSsn_clim_maxLag_%d_%s_%d_%d.mat', lag, prod, year_st_pheno, year_ed_pheno);
        save(fullfile(dir_sens, file_preSsn), '-fromstruct', s_preSsn, '-v7.3');
    end
    elapsed_time = toc(start_time);
    fprintf('[Info] Processing time: %.2f min\n', elapsed_time / 60);

end

