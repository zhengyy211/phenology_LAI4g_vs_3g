%% ========================================================================
%  File:    code_2_2_prep_vi.m
%  Author:  Yaoyao Zheng
%  Email:   zhengyaoyao@stu.pku.edu.cn
%  MATLAB:  R2026a
%
%  Description:
%  Preprocess VI products for subsequent analysis. Tasks include:
%    1) Reading the original VI product, Cropping to the study region (≥30°N).
%    2) Filling missing values with zero.
%    3) Detecting and correcting spikes during the thermal non-growing season.
%    4) Saving the processed dataset.
% =========================================================================

clear;

%% ------ STEP 1: Load configuration ------
run('code_0_0_config.m');

% Define directories for raw VI data and processed outputs
dir_raw   = cfg.path.raw;
dir_prep  = cfg.path.prep;

% List of VI products to process
products  = [cfg.vi.products, 'modis_lai'];

% Spatial range: Northern Hemisphere
lat_limit = cfg.space.lat_limit;
n_rows    = cfg.space.n_rows;     % Target rows after resampling
n_cols    = cfg.space.n_cols;     % Target columns after resampling

max_workers = min(16, feature('numcores'));
pool = gcp('nocreate');
if isempty(pool)
    parpool(max_workers);
elseif pool.NumWorkers ~= max_workers
    delete(pool);
    parpool(max_workers);
end

%% ------ STEP 2: Process VI products ------
for prod = products

    % prod = products(7);
    fprintf('\n[Info] ====== Processing %s ======\n', prod);
    start_time = tic;

    dir_raw_vi = fullfile(dir_raw, prod);
    year_start = cfg.vi.(prod).year_start;
    year_end   = cfg.vi.(prod).year_end;
    year_list  = year_start:1:year_end;
    doy_list   = cfg.vi.(prod).doy_list;
    n_years    = numel(year_list);
    n_doys     = numel(doy_list);

    % Determine latitude range
    res      = cfg.vi.(prod).spatial_res;
    lat_vals = ((90 - res/2):-res:(-90 + res/2));
    row_st   = find(lat_vals <= lat_limit(2), 1, 'first');
    row_ed   = find(lat_vals >= lat_limit(1), 1,  'last');

    %% ------ STEP 2.1: Read VI and QC information ------
    % ------ GIMMS ------
    if contains(prod, 'gimms') | strcmpi(prod, 'mc_ndvi')

        data_vi      = cell(n_years, 1);
        data_weights = cell(n_years, 1);

        parfor idx_yr = 1:n_years

            yr = year_list(idx_yr);
            fprintf('[Info] Year %d...\n', yr);

            % Read data
            switch lower(prod)
                case 'gimms_lai4g'
                    [vi, weights] = read_gimms_lai4g(dir_raw_vi, yr);
                case 'gimms_lai3g'
                    [vi, weights] = read_gimms_lai3g(dir_raw_vi, yr);
                case 'gimms_ndvi4g'
                    [vi, weights] = read_gimms_ndvi4g(dir_raw_vi, yr);
                case 'gimms_ndvi3g'
                    [vi, weights] = read_gimms_ndvi3g(dir_raw_vi, yr);
                case 'mc_ndvi'
                    [vi, weights] = read_mc_ndvi(dir_raw_vi, yr);
            end

            % Crop to the Northern Hemisphere region and reshape to column vector
            vi      =      vi(row_st:row_ed, :, :);
            weights = weights(row_st:row_ed, :, :);

            % Perform Spatial Interpolation
            [r, c, t] = size(vi, [1, 2, 3]);
            if r == n_rows && c == n_cols
                % do nothing
            elseif r < n_rows && c < n_cols
                cell_vi = cell(t, 1);
                cell_w  = cell(t, 1);
                for i_t = 1:t
                    cell_vi{i_t} = imresize(vi(:, :, i_t), [n_rows, n_cols], 'bilinear');
                    cell_w{i_t}  = imresize(vi(:, :, i_t), [n_rows, n_cols], 'nearest');
                end
                vi      = cat(3, cell_vi{:});
                weights = cat(3, cell_w{:});
            elseif r > n_rows && c > n_cols
                vi      = imresize_block_mean(vi,      n_rows, n_cols);
                weights = imresize_block_mean(weights, n_rows, n_cols);
            else
                error('Check data size mismatch.');
            end

            data_vi{idx_yr}      = reshape(vi,      n_rows * n_cols, t);
            data_weights{idx_yr} = reshape(weights, n_rows * n_cols, t);

        end

        % Concatenate data for all time steps into a [pixel x time] matrix
        data_vi      = horzcat(data_vi{:});
        data_weights = horzcat(data_weights{:});

        clearvars vi weights

    else

        data_vi = cell(n_years, 1);

        parfor idx_yr = 1:n_years

            yr = year_list(idx_yr);
            fprintf('[Info] Year %d...\n', yr);

            % Read data
            switch lower(prod)
                case 'glass_lai'
                    vi = read_glass_lai(dir_raw_vi, yr);
                case 'yuan_lai'
                    vi = read_yuan_lai(dir_raw_vi, yr);
                case 'modis_lai'
                    vi = read_mod15a2h(dir_raw_vi, yr);
            end

            % Crop to the Northern Hemisphere region and reshape to column vector
            vi = vi(row_st:row_ed, :, :);

            % Perform Spatial Interpolation
            [r, c, t] = size(vi, [1, 2, 3]);
            if r == n_rows && c == n_cols
                % do nothing
            elseif r < n_rows && c < n_cols
                cell_vi = cell(t, 1);
                for i_t = 1:t
                    cell_vi{i_t} = imresize(vi(:, :, i_t), [n_rows, n_cols], 'bilinear');
                end
                vi = cat(3, cell_vi{:});
            elseif r > n_rows && c > n_cols
                vi = imresize_block_mean(vi, n_rows, n_cols);
            else
                error('Check data size mismatch.');
            end

            data_vi{idx_yr} = reshape(vi, n_rows * n_cols, size(vi, 3));

        end

        % Concatenate data for all time steps into a [pixel x time] matrix
        data_vi = horzcat(data_vi{:});

        clearvars vi

    end


    %% ------ STEP 2.2: Apply vegetation mask ------
    % Calculate the multi-year average VI for each pixel and retain pixels with
    % an average VI > 0.1
    data_reshape = reshape(data_vi, size(data_vi, 1), n_doys, n_years);
    yearly_max   = squeeze(max(data_reshape, [], 2, 'omitmissing'));
    yearly_min   = squeeze(min(data_reshape, [], 2, 'omitmissing'));
    mask_ampl    = all((yearly_max - yearly_min) > 0.1, 2);
    clearvars data_reshape yearly_max yearly_min

    mask_mean = sum(data_vi, 2, 'omitmissing') / size(data_vi, 2) > 0.1;
    mask_vi   = mask_mean & mask_ampl;
    data_vi   = data_vi(mask_vi, :);

    if contains(prod, 'gimms')
        data_weights = data_weights(mask_vi, :);
    else
        data_weights = ones(size(data_vi), 'single') * 0.2;
        data_weights(~isnan(data_vi)) = 1;
    end

    if contains(prod, 'gimms_lai') | strcmpi(prod, 'modis_lai')
        file_vi_raw  = sprintf('data_raw_%s_%d_%d.mat', prod, year_start, year_end);
        save(fullfile(dir_raw, file_vi_raw), 'data_vi', 'mask_vi', '-v7.3');
    end

    %% ------ STEP 2.3: Fill missing values ------
    mask_missing = isnan(data_vi);
    data_vi(mask_missing) = 0;

    %% ------ STEP 2.4: Spike detection and correction ------
    % Load thermal growing season
    % thermal_gs : [all pixel x time] logical mask for temperature > 0°C (thermal growing season)
    file_temp = sprintf('data_thermal_gs_%s_%d_%d.mat', prod, year_start, year_end);
    load(fullfile(dir_prep, file_temp), 'thermal_gs');
    thermal_gs = thermal_gs(mask_vi, :);

    % A spike is defined as when the absolute difference between the observed VI
    % and the 3-point moving average exceeds 10% of the pixel's maximum value.
    [data_vi, mask_spikes] = remove_spikes(data_vi, data_weights, thermal_gs);
    clearvars thermal_gs


    %% ------ STEP 2.5: Update quality weights ------
    % Further downweight filled or spike-corrected pixels
    data_weights(mask_missing) = 0.2;
    data_weights(mask_spikes)  = 0.2;


    %% ------ STEP 2.6: Save processed data ------
    % mask_vi      : [all pixel   x    1] Logical mask for valid vegetation pixels
    % data_vi      : [valid pixel x time] VI values
    % data_weights : [valid pixel x time] Quality weights corresponding to data_vi
    % n_pixels     : Number of valid pixels
    n_pixels = sum(mask_vi);
    file_vi  = sprintf('data_%s_%d_%d.mat', prod, year_start, year_end);
    save(fullfile(dir_prep, file_vi), 'data_vi', 'data_weights', 'mask_vi', 'n_pixels', '-v7.3');

    %% ------ STEP 2.7: Plot regional mean VI for visual check ------
    mean_vi = mean(data_vi, 1, 'omitmissing');

    figure('Name', sprintf('[%s] Mean VI', prod), 'Color', 'w');
    plot(mean_vi, '-b', 'LineWidth', 1.5);
    ylabel('VI');
    title(sprintf('[%s] Regional mean VI', prod));
    grid on;

    clearvars data_vi data_weights mask_vi

    elapsed_time = toc(start_time);
    fprintf('[Info] Processing time for %s: %.2f min\n', prod, elapsed_time / 60);

end

delete(gcp('nocreate'));