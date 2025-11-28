%% ========================================================================
%  File:        code_2_5_prep_mcd12q2.m
%  Author:      Yaoyao Zheng
%  Email:       zhengyaoyao@stu.pku.edu.cn
%  MATLAB:      R2025b
%  Description: This script prepares MODIS MCD12Q2 phenology data for
%               subsequent analysis of vegetation phenology.
%               1) Reading MODIS TIFF files for each year and phase.
%               2) Converting raw time units ('days since 1970') to
%                  day-of-year (DOY).
%               3) Aggregating spatial resolution via block circ mean to
%                  match the study grid.
%               4) Ensuring phenological events do not cross year boundaries 
%                  by applying seasonal shifts to the DOY values.
% =========================================================================

clear;

%% ------ STEP 1: Load configuration ------
run('code_0_0_config.m');

% Define directories for raw MODIS data and processed outputs
dir_modis  = fullfile(cfg.path.raw, 'mcd12q2');
dir_prep   = cfg.path.prep;
dir_pheno  = cfg.path.pheno;

% Define study period for MODIS
year_start = cfg.modis.year_start;
year_end   = cfg.modis.year_end;
year_list  = year_start:year_end;
n_years    = numel(year_list);

% Spatial grid parameters
n_rows     = cfg.space.n_rows;
n_cols     = cfg.space.n_cols;

% Phenological phases
phases = cfg.pheno.modis.phases;

%% ------ STEP 2: Process phenology data ------
s_pheno = struct;

mask_pheno = true(n_rows * n_cols, 1);

for phase = phases

    fprintf('\n[Info] ====== Processing %s ======\n', phase);

    cell_pheno = cell(n_years, 1);

    for idx_yr = 1:n_years

        current_year = year_list(idx_yr);
        fprintf('[Info] Year %d...\n', current_year);

        %% ------ STEP 2.1: Read East/West tiles ------
        file_pattern = sprintf('MCD12Q2_%s_1_%04d', phase, current_year);
        data_W = imread(fullfile(dir_modis, [file_pattern, '_W.tif']));  % Western tile
        data_E = imread(fullfile(dir_modis, [file_pattern, '_E.tif']));  % Eastern tile

        %% ------ STEP 2.2: Convert to Day-of-Year (DOY) ------
        % convert from 'days_since_1970' to DOY
        data_doy = [get_mcd12q2_doy(data_W), get_mcd12q2_doy(data_E)];
        clearvars data_W data_E

        %% ------ STEP 2.3: Aggregate to target spatial resolution ------
        data_resized = imresize_block_circ_mean(data_doy, n_rows, n_cols);
        clearvars data_doy

        cell_pheno{idx_yr} = reshape(data_resized, [], 1);

    end
    s_pheno.(phase) = horzcat(cell_pheno{:}); % [pixel × time]

    mask_pheno = mask_pheno & any(~isnan(s_pheno.(phase)), 2);

end

for phase = phases
    s_pheno.(phase) = s_pheno.(phase)(mask_pheno, :);
end
s_pheno.mask_pheno = mask_pheno;

%% ------ STEP 3: Save processed phenology data ------
file_mcd12q2 = sprintf('data_mcd12q2_%d_%d.mat', year_start, year_end);
save(fullfile(dir_prep, file_mcd12q2), '-fromstruct', s_pheno, '-v7.3')


%% ------ STEP 4: Apply seasonal shifts to avoid year-crossing ------
% In this step, we apply seasonal shifts to ensure phenology events (like Greenup)
% do not span across different years.

file_mcd12q2 = sprintf('data_mcd12q2_%04d_%04d.mat', year_start, year_end);
s_pheno = load(fullfile(dir_prep, file_mcd12q2));

% Compute the average peak day-of-year (DOY) for each pixel across years
peak_doy      = s_pheno.Peak;
peak_doy_avg  = cellfun(@(x) circ_mean_doy(x), num2cell(peak_doy, 2), 'UniformOutput', false);
peak_doy_avg  = vertcat(peak_doy_avg{:});

% DOY ranges for different seasons
doys_spring   = 60:151;           % ~March–May
doys_summer   = 152:243;          % ~June–August
doys_autumn   = 244:334;          % ~September–November
doys_winter   = [335:365, 1:59];  % ~December–February

% Determine vegetation type by peak DOY
veg_type = NaN(size(peak_doy_avg));
veg_type(ismember(peak_doy_avg, doys_spring)) = 1;
veg_type(ismember(peak_doy_avg, doys_summer)) = 2;
veg_type(ismember(peak_doy_avg, doys_autumn)) = 3;
veg_type(ismember(peak_doy_avg, doys_winter)) = 4;

% Initialize shifted struct
s_pheno_shifted            = struct;
s_pheno_shifted.veg_type   = veg_type;
s_pheno_shifted.mask_pheno = s_pheno.mask_pheno;

% Define seasonal shift in months for each vegetation type to avoid year-crossing
season_shift_months = struct(...
    'spring',  3, ...   % spring -> shift +3 months
    'summer',  0, ...   % Summer -> no shift
    'autumn', -3, ...   % autumn -> shift -3 months
    'winter',  6 ...    % winter -> shift +6 months
    );

% Apply the seasonal shift to the phenology data (DOY)
for phase = phases
    data_phase = s_pheno.(phase);  % [n_pixel x n_year]
    data_shifted = apply_seasonal_shift_to_doy( ...
                   data_phase, veg_type, season_shift_months);
    s_pheno_shifted.(phase) = data_shifted;
end

% Save the shifted phenology data
file_mcd12q2 = sprintf('data_pheno_raw_mcd12q2_%d_%d.mat', year_start, year_end);
save(fullfile(dir_pheno, file_mcd12q2), '-fromstruct', s_pheno_shifted, '-v7.3')