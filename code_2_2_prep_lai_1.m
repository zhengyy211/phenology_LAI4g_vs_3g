%% ========================================================================
%  File:        code_2_2_prep_lai_1.m
%  Author:      Yaoyao Zheng
%  Email:       zhengyaoyao@stu.pku.edu.cn
%  MATLAB:      R2025b
%  Description: Preprocess LAI products 
%               GIMMS LAI4g, GIMMS LAI3g, GLASS LAI, MODIS LAI (MOD15A2H)
%               1) Load data and convert to 1/12° resolution if necessary
%               2) Crop to Northern Hemisphere (≥30°N) and apply vegetation mask
%               3) Save preprocessed LAI and QC data
% =========================================================================

clear;

%% ------ STEP 1: Load configuration ------
run('code_0_0_config.m');

% Define directories for raw LAI data and processed outputs
dir_raw    = cfg.path.raw;
dir_prep   = cfg.path.prep;

% Time range
year_start = cfg.lai.year_start;
year_end   = cfg.lai.year_end; 
year_list  = year_start:1:year_end;
doy_list   = cfg.lai.glass_yuan.doy_list;

% Spatial range: Northern Hemisphere
row_st     = cfg.space.row_st;     % Northern Hemisphere start row
row_ed     = cfg.space.row_ed;     % Northern Hemisphere end row
lat_limit  = cfg.space.lat_limit;  
n_rows     = cfg.space.n_rows;     % Target rows after resampling
n_cols     = cfg.space.n_cols;     % Target columns after resampling

%% ------ STEP 2: GIMMS LAI (LAI3g / LAI4g) ------
% Generate half-monthly combinations
months = 1:12;
halves = [1 2]; % 1 = a, 2 = b
[yy, mm, hh] = ndgrid(year_list, months, halves);
combs_time   = table(yy(:), mm(:), hh(:), 'VariableNames', {'year', 'month', 'half'});
combs_time   = sortrows(combs_time, {'year', 'month', 'half'});
clearvars months halves yy mm hh

for lai_prod = ["gimms_lai4g", "gimms_lai3g"]

    fprintf('\n[Info] ====== Reading %s ======\n', lai_prod);
    start_time = tic;

    %% ------ STEP 2.1: Read LAI and QC data for each half-month ------
    % For each half-month time step, read LAI and QC data, crop to the Northern 
    % Hemisphere region, and store in a cell array, which is later concatenated 
    % into a matrix of size [pixel x time].

    data_lai = cell(height(combs_time), 1);
    data_qc  = cell(height(combs_time), 1);

    parfor idx_comb = 1:height(combs_time)

        current_year  = combs_time.year(idx_comb);
        current_month = combs_time.month(idx_comb);
        current_half  = combs_time.half(idx_comb);

        if current_month == 1 && current_half == 1
            fprintf('[Info] Year %d...\n', current_year);
        end

        % Read LAI and QC data
        if lai_prod == "gimms_lai3g"
            [lai, qc] = read_gimms_lai3g(fullfile(dir_raw, lai_prod), current_year, current_month, current_half);
        elseif lai_prod == "gimms_lai4g"
            [lai, qc] = read_gimms_lai4g(fullfile(dir_raw, lai_prod), current_year, current_month, current_half);
        end

        % Plot first image for visual check
        if current_year == year_list(1) && current_month == 1 && current_half == 1
            figure('Name', sprintf('[%s] LAI', lai_prod), 'Color', 'w');
            h = imagesc(lai);
            cb = colorbar('southoutside');
            cb.Label.String = 'LAI (m^2/m^2)';
            cb.FontWeight = 'bold';
            h.AlphaData = ~isnan(lai);
            title(sprintf('[%s] LAI', lai_prod));

            figure('Name', sprintf('[%s] QC', lai_prod), 'Color', 'w');
            h = imagesc(qc);
            cb = colorbar('southoutside');
            cb.Label.String = 'QC';
            cb.FontWeight = 'bold';
            h.AlphaData = ~isnan(qc);
            title(sprintf('[%s] QC', lai_prod));
        end

        % Crop to the Northern Hemisphere region and reshape to column vector
        data_lai{idx_comb} = reshape(lai(row_st:row_ed, :), [], 1);
        data_qc{idx_comb}  = reshape( qc(row_st:row_ed, :), [], 1);

    end

    % Concatenate data for all time steps into a [pixel x time] matrix
    data_lai = horzcat(data_lai{:});
    data_qc  = horzcat(data_qc{:});

    %% ------ STEP 2.2: Apply vegetation mask ------
    % Calculate the multi-year average LAI for each pixel and retain pixels with 
    % an average LAI > 0.1
    mask_lai = sum(data_lai, 2, 'omitmissing') / size(data_lai, 2) > 0.1;
    data_lai = data_lai(mask_lai, :);
    data_qc  = data_qc(mask_lai, :);

    %% ------ STEP 2.3: Save raw data ------
    % Save raw data, including:
    % mask_lai     : [all pixel   x    1] logical mask for valid vegetation pixels
    % data_lai     : [valid pixel x time] filled and corrected LAI values
    % data_qc      : [valid pixel x time] QC values
    file_raw = sprintf('data_raw_%s_%d_%d.mat', lai_prod, year_start, year_end);
    save(fullfile(dir_raw, file_raw), 'data_lai', 'data_qc', 'mask_lai', '-v7.3');

    clearvars data_lai data_qc mask_lai lai qc

    elapsed_time = toc(start_time);
    fprintf('[Info] Processing time for %s: %.2f min\n', lai_prod, elapsed_time / 60);

end


%% ------ STEP 3: GLASS LAI ------
% Generate combinations of year and doy
[yy, dd]   = ndgrid(year_list, doy_list);
combs_time = table(yy(:), dd(:), 'VariableNames', {'year', 'doy'});
combs_time = sortrows(combs_time, {'year', 'doy'});
clearvars doy_list yy dd

% Row indices for latitude range
glass_spatial_res = 0.05;  % GLASS dataset has 0.05° resolution
glass_lat_vals    = ((90 - glass_spatial_res/2):-glass_spatial_res:(-90 + glass_spatial_res/2));
glass_row_st      = find(glass_lat_vals <= lat_limit(2), 1, 'first');
glass_row_ed      = find(glass_lat_vals >= lat_limit(1), 1,  'last');

fprintf('\n[Info] ====== Reading GLASS LAI ======\n');
start_time = tic;

%% ------ STEP 3.1: Read LAI data for each half-month ------

data_lai = cell(height(combs_time), 1);

parfor idx_comb = 1:height(combs_time)

    current_year = combs_time.year(idx_comb);
    current_doy  = combs_time.doy(idx_comb);

    if current_doy == 1
        fprintf('[Info] Year %d...\n', current_year);
    end

    % Read LAI data
    lai = read_glass_lai(fullfile(dir_raw, 'glass_lai'), current_year, current_doy);

    % Plot first image for visual check
    if current_year == year_list(1) && current_doy == 1
        figure('Name', '[GLASS] LAI', 'Color', 'w');
        h = imagesc(lai);
        cb = colorbar('southoutside');
        cb.Label.String = 'LAI (m^2/m^2)';
        cb.FontWeight = 'bold';
        h.AlphaData = ~isnan(lai);
        title('[GLASS] LAI');
    end

    % Crop to the Northern Hemisphere region
    lai = lai(glass_row_st:glass_row_ed, :);
    lai = imresize_block_mean(lai, n_rows, n_cols);  % Convert to target resolution
    data_lai{idx_comb} = reshape(lai, [], 1);

end

% Concatenate data for all time steps into a [pixel x time] matrix
data_lai = horzcat(data_lai{:});

%% ------ STEP 3.2: Apply vegetation mask ------
% Calculate the multi-year average LAI for each pixel and retain pixels with
% an average LAI > 0.1
mask_lai = sum(data_lai, 2, 'omitmissing') / size(data_lai, 2) > 0.1;
data_lai = data_lai(mask_lai, :);

% Define QC: 0 = good, 3 = missing
data_qc = zeros(size(data_lai));
data_qc(isnan(data_lai)) = 3;

%% ------ STEP 3.3: Save raw data ------
% Save raw data, including:
% mask_lai     : [all pixel   x    1] logical mask for valid vegetation pixels
% data_lai     : [valid pixel x time] filled and corrected LAI values
% data_qc      : [valid pixel x time] QC values
file_raw = sprintf('data_raw_glass_lai_%d_%d.mat', year_start, year_end);
save(fullfile(dir_raw, file_raw), 'data_lai', 'data_qc', 'mask_lai', '-v7.3');

clearvars data_lai data_qc mask_lai lai qc

elapsed_time = toc(start_time);
fprintf('[Info] Processing time for GLASS LAI: %.2f min\n', elapsed_time / 60);


%% ------ STEP 4: Yuan LAI v6.1 ------
fprintf('\n[Info] ====== Reading Yuan LAI ======\n');
start_time = tic;

year_start = cfg.modis.year_start;
year_end   = cfg.modis.year_end; 
year_list  = year_start:1:year_end;
n_years    = numel(year_list);

% Row indices for latitude range
yuan_spatial_res = 0.1;  % Yuan dataset has 0.1° resolution
yuan_lat_vals    = ((90 - yuan_spatial_res/2):-yuan_spatial_res:(-90 + yuan_spatial_res/2));
yuan_row_st      = find(yuan_lat_vals <= lat_limit(2), 1, 'first');
yuan_row_ed      = find(yuan_lat_vals >= lat_limit(1), 1,  'last');

%% ------ STEP 4.1: Read LAI data for each year ------
data_lai = cell(n_years, 1);
parfor idx_year = 1:n_years

    current_year = year_list(idx_year);
    fprintf('[Info] Year %d...\n', current_year);

    % Read LAI data
    lai = read_yuan_lai(fullfile(dir_raw, 'yuan_lai'), current_year);

    % Plot first image for visual check
    if current_year == year_list(1)
        figure('Name', '[Yuan] LAI', 'Color', 'w');
        h = imagesc(lai(:,:,1));
        cb = colorbar('southoutside');
        cb.Label.String = 'LAI (m^2/m^2)';
        cb.FontWeight = 'bold';
        h.AlphaData = ~isnan(lai(:,:,1));
        title('[Yuan] LAI');
    end

    lai = lai(yuan_row_st:yuan_row_ed, :, :);
    lai_resized = cell(size(lai, 3), 1);
    for ii = 1:size(lai, 3)
        lai_resized{ii} = imresize(lai(:, :, ii), [n_rows, n_cols], 'bilinear');
    end
    lai_resized = cat(3, lai_resized{:});
    
    data_lai{idx_year} = reshape(lai_resized, n_rows * n_cols, 46);

end

% Concatenate data for all time steps into a [pixel x time] matrix
data_lai = horzcat(data_lai{:});

%% ------ STEP 4.2: Apply vegetation mask ------
% Calculate the multi-year average LAI for each pixel and retain pixels with
% an average LAI > 0.1
mask_lai = sum(data_lai, 2, 'omitmissing') / size(data_lai, 2) > 0.1;
data_lai = data_lai(mask_lai, :);

% Define QC: 0 = good, 3 = missing
data_qc = zeros(size(data_lai));
data_qc(isnan(data_lai)) = 3;

%% ------ STEP 4.3: Save raw data ------
% Save raw data, including:
% mask_lai     : [all pixel   x    1] logical mask for valid vegetation pixels
% data_lai     : [valid pixel x time] filled and corrected LAI values
% data_qc      : [valid pixel x time] QC values
file_raw = sprintf('data_raw_yuan_lai_%d_%d.mat', year_start, year_end);
save(fullfile(dir_raw, file_raw), 'data_lai', 'data_qc', 'mask_lai', '-v7.3');

clearvars data_lai data_qc mask_lai lai qc

elapsed_time = toc(start_time);
fprintf('[Info] Processing time for Yuan LAI: %.2f min\n', elapsed_time / 60);


%% ------ STEP 5: MOD15A2H ------

fprintf('\n[Info] ====== Reading MODIS LAI (MOD15A2H) ======\n');
start_time = tic;

year_start = cfg.modis.year_start;
year_end   = cfg.modis.year_end; 
year_list  = year_start:1:year_end;
n_years    = numel(year_list);

%% ------ STEP 5.1: Read LAI data for each year ------
data_lai = cell(n_years, 1);
parfor idx_year = 1:n_years

    current_year = year_list(idx_year);
    fprintf('[Info] Year %d...\n', current_year);

    % Read LAI data
    lai = read_mod15a2h(fullfile(dir_raw, 'mod15a2h'), current_year);

    % Plot first image for visual check
    if current_year == year_list(1)
        figure('Name', '[MOD15A2H] LAI', 'Color', 'w');
        h = imagesc(lai(:,:,1));
        cb = colorbar('southoutside');
        cb.Label.String = 'LAI (m^2/m^2)';
        cb.FontWeight = 'bold';
        h.AlphaData = ~isnan(lai(:,:,1));
        title('[MOD15A2H] LAI');
    end

    data_lai{idx_year} = reshape(lai, size(lai, 1) * size(lai, 2), 46);

end

% Concatenate data for all time steps into a [pixel x time] matrix
data_lai = horzcat(data_lai{:});

%% ------ STEP 5.2: Apply vegetation mask ------
% Calculate the multi-year average LAI for each pixel and retain pixels with
% an average LAI > 0.1
mask_lai = sum(data_lai, 2, 'omitmissing') / size(data_lai, 2) > 0.1;
data_lai = data_lai(mask_lai, :);

% Define QC: 0 = good, 3 = missing
data_qc = zeros(size(data_lai));
data_qc(isnan(data_lai)) = 3;

%% ------ STEP 5.3: Save raw data ------
% Save raw data, including:
% mask_lai     : [all pixel   x    1] logical mask for valid vegetation pixels
% data_lai     : [valid pixel x time] filled and corrected LAI values
% data_qc      : [valid pixel x time] QC values
file_raw = sprintf('data_raw_modis_lai_%d_%d.mat', year_start, year_end);
save(fullfile(dir_raw, file_raw), 'data_lai', 'data_qc', 'mask_lai', '-v7.3');

clearvars data_lai data_qc mask_lai lai qc

elapsed_time = toc(start_time);
fprintf('[Info] Processing time for MODIS LAI (MOD15A2H): %.2f min\n', elapsed_time / 60);