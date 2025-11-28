%% ========================================================================
%  File:        code_4_2_sens_coef.m
%  Author:      Yaoyao Zheng
%  Email:       zhengyaoyao@stu.pku.edu.cn
%  MATLAB:      R2025b
%  Description: Computes pixel-wise sensitivity of phenological transition dates
%               to pre-season climate variables.
%               For each LAI product and phenological phase:
%               1) Extracts cleaned phenology data and corresponding climate data,
%               2) Computes climate time series at multiple pre-season lags,
%               3) Identifies the optimal lag maximizing correlation with phenology,
%               4) Stores both the optimal lag and aligned pre-season climate data.
%               Outputs are saved for subsequent climate sensitivity calculations.
% =========================================================================

clear;

%% ------ STEP 1: Load configuration ------
run('code_0_0_config.m');

% Define directories for sensitivity results
dir_sens     = cfg.path.sens;

% Define variables for processing
lai_products = cfg.lai.products;
products     = [lai_products, "mcd12q2"];

% Climate variables to analyze
vars_clim    = cfg.sens.vars_clim;
n_clims      = numel(vars_clim);

max_lag_mths = cfg.sens.max_lag_mths;  % Maximum number of lag periods to consider

% Minimum valid observations required for the regression
% min_obs      = n_clims + 2;
min_obs      = 10;

% Spatial / grid info
n_rows       = cfg.space.n_rows;
n_cols       = cfg.space.n_cols;
R            = cfg.space.R;

%% ------ STEP 2: Loop through each product (LAI or MODIS) ------
for prod = products

    fprintf('\n[Info] ====== Processing %s ======\n', prod);
    start_time = tic;

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

    for lag = max_lag_mths
        file_suffix = sprintf('maxLag_%d_%s_%d_%d', lag, prod, year_start, year_end);
        file_preSsn = sprintf('data_preSsn_clim_%s.mat', file_suffix);
        s_preSsn    = load(fullfile(dir_sens, file_preSsn));

        s_sens           = struct;
        s_sens.mask_sens = s_preSsn.mask_pheno;
        n_pixels         = sum(s_sens.mask_sens);
        mask_1d          = s_sens.mask_sens;

        for phase = phases

            fprintf('[Info] %s ...\n', phase);

            %% ------ STEP 2.1: Load the corresponding pre-season climate data ------
            data_clim = cell(n_clims, 1);
            for idx_clim = 1:n_clims
                data_clim{idx_clim} = s_preSsn.(phase).(sprintf('preSsn_%s', vars_clim(idx_clim)));
            end
            data_clim = cat(3, data_clim{:});

            % Identify valid climate data (not NaN) for each pixel and time
            % Check if climate data is valid at any time
            valid_clim = squeeze(any(~isnan(data_clim), 2));
            % Ensure data is valid for all climate variables
            valid_clim = all(valid_clim, 2);
            data_clim  = num2cell(data_clim, [2 3]);
            data_clim  = cellfun(@(x) squeeze(x), data_clim, 'UniformOutput', false);

            %% ------ STEP 2.2: Load phenology data ------
            data_pheno  = s_preSsn.(phase).pheno;
            valid_pheno = sum(~isnan(data_pheno), 2) >= min_obs;
            data_pheno  = num2cell(s_preSsn.(phase).pheno, 2);

            % Apply masks for valid pixels (both phenology and climate)
            valid_pixels   = valid_pheno & valid_clim;
            n_valid_pixels = sum(valid_pixels);
            data_pheno     = data_pheno(valid_pixels);
            data_clim      = data_clim(valid_pixels);

            %% ------ STEP 2.3: Pre-season climate sensitivity ------
            % Initialize containers for regression results
            cell_beta    = cell(n_valid_pixels, 1);
            cell_p_beta  = cell(n_valid_pixels, 1);
            cell_R2      = cell(n_valid_pixels, 1);
            parfor pixel = 1:n_valid_pixels
                [cell_beta{pixel}, cell_p_beta{pixel}, cell_R2{pixel}] = ...
                    calc_clim_sens(data_pheno{pixel}, data_clim{pixel}, min_obs);
            end

            %% ------ STEP 2.4: Organize the results ------
            beta_full = NaN(n_pixels, n_clims, 'single');
            pval_full = NaN(n_pixels, n_clims, 'single');
            R2_full   = NaN(n_pixels, 1, 'single');

            % Assign the calculated values for valid pixels
            beta_full(valid_pixels, :) = vertcat(cell_beta{:});
            pval_full(valid_pixels, :) = vertcat(cell_p_beta{:});
            R2_full(valid_pixels)      = vertcat(cell_R2{:});

            % Store the results in the sensitivity structure
            s_sens.(phase).beta = beta_full;
            s_sens.(phase).pval = pval_full;
            s_sens.(phase).R2   = R2_full;


            %% ------ STEP 2.5: Write results to GeoTIFFs ------
            for idx_clim = 1:n_clims
                file_beta = fullfile(dir_sens, sprintf('map_sens_beta_%s_%s_%s.tif', phase, vars_clim(idx_clim), file_suffix));
                write_map_to_tiff(beta_full(:, idx_clim), mask_1d, n_rows, n_cols, R, file_beta);

                file_pval = fullfile(dir_sens, sprintf('map_sens_pval_%s_%s_%s.tif', phase, vars_clim(idx_clim), file_suffix));
                write_map_to_tiff(pval_full(:, idx_clim), mask_1d, n_rows, n_cols, R, file_pval);
            end
        end

        %% ------ STEP 2.6: Save the final sensitivity results -----------------
        file_sens = sprintf('data_sens_%s.mat', file_suffix);
        save(fullfile(dir_sens, file_sens), '-fromstruct', s_sens, '-v7.3')
    end

    elapsed_time = toc(start_time);
    fprintf('[Info] Processing time: %.2f min\n',elapsed_time / 60);

end