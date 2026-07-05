%% ========================================================================
%  File:    code_2_4_prep_phenocam.m
%  Author:  Yaoyao Zheng
%  Email:   zhengyaoyao@stu.pku.edu.cn
%  MATLAB:  R2026a
%
%  Description:
%  This script processes and prepares the Phenocam data for phenological 
%  analysis. It includes steps for loading metadata, extracting phenology 
%  phases, cleaning, and spatial gridding.
% =========================================================================

clear;

%% ------ STEP 1: Load configuration ------
run('code_0_0_config.m');

% Define directories for raw phenocam data and processed outputs
dir_raw      = fullfile(cfg.path.raw, 'PhenoCam');
dir_prep     = cfg.path.prep;

% Time range for processing
year_start   = cfg.time.year_start;
year_end     = cfg.time.year_end;
t_start      = datetime(year_start,  1,  1);
t_end        = datetime(year_end,   12, 31);

phase_lvls   = cfg.pheno.phase_lvls;
phase_lbls   = cfg.pheno.phase_lbls;

% Spatial / grid info
spatial_res  = cfg.space.spatial_res;
lat_min      = cfg.space.lat_limit(1);

% Set Z-score threshold for outlier detection
z_score_thr = cfg.pheno.z_score_thr;


%% ------ STEP 2: Unzip dataset ------
% Unzip all data files if they are in compressed (.zip) format
file_list = dir(fullfile(dir_raw, '*.zip'));
for i_list = 1:size(file_list, 1)
    file_zip = fullfile(file_list(i_list).folder, file_list(i_list).name);

    [~, name_zip, ~] = fileparts(file_list(i_list).name);
    dir_unzip = fullfile(dir_raw, name_zip);

    if ~exist(dir_unzip, 'dir')
        fprintf('[Info] Unzipping: %s\n', file_zip);
        unzip(file_zip, dir_unzip);
    end
end


%% ------ STEP 3: Load metadata ------ 
% Read the metadata for each Phenocam dataset folder
folders = dir(fullfile(dir_raw, 'PhenoCam_V3_*'));
folders = folders([folders.isdir]);

n_folders = numel(folders);

cell_meta = cell(n_folders, 1);

% Load site and ROI metadata for each folder
for i_folder = 1:n_folders

    folder = fullfile(dir_raw, folders(i_folder).name);

    % Read site metadata
    site_meta_file = dir(fullfile(folder, 'data_record_1', '*_meta.txt'));
    site_meta_path = fullfile(site_meta_file.folder, site_meta_file.name);
    site_meta      = read_phenocam_site_meta(site_meta_path);

    % Read ROI metadata
    roi_meta_file  = dir(fullfile(folder, 'data_record_2', '*_roi.csv'));
    roi_meta_path  = fullfile(roi_meta_file.folder, roi_meta_file.name);
    roi_meta       = read_phenocam_roi_meta(roi_meta_path);
    
    % Add folder path for reference
    roi_meta.folder = string(folder);

    % Combine metadata
    cell_meta{i_folder} = [site_meta, roi_meta];

end
% Concatenate all metadata
data_meta = vertcat(cell_meta{:});


%% ------ STEP 4: Filter and process metadata ------ 
% Filter metadata based on latitude
data_meta = data_meta(data_meta.lat >= lat_min, :);

% Keep only the sites that have data within the specified time range
overlap_start = max(data_meta.date_start, t_start);
overlap_end   = min(data_meta.date_end,   t_end);
overlap_years = years(overlap_end - overlap_start);
data_meta     = data_meta(overlap_years >= 0, :);

data_meta     = data_meta(data_meta.primary_veg_type == data_meta.roi_veg_type | ...
                          data_meta.primary_veg_type == "", :);

%% ------ STEP 5: Extract the phenology data ------
n_sites = height(data_meta);
cell_pheno = cell(n_sites, 1);

parfor i = 1:n_sites
    % Directory for the simplified daily data for each site
    folder    = fullfile(data_meta.folder(i), 'data_record_7');
    csv_files = dir(fullfile(folder, '*_simplified_1day.csv'));

    if isempty(csv_files)
        warning('No simplified_1day CSV found for site %s roi %d', ...
            data_meta.sitename(i), data_meta.roi_id(i));
        continue;
    end

    % Read the CSV file containing the simplified daily data
    csv_path  = fullfile(csv_files(1).folder, csv_files(1).name);
    data_1day = readtable(csv_path);
    
    % Extract phenology phases for the current site
    site_pheno = extract_pheno_phenocam(data_1day, cfg);
    
    % Reshape the phenology data into a long format (melt the data)
    id_col     = 'year';
    value_vars = site_pheno.Properties.VariableNames;
    value_vars = setdiff(value_vars, id_col);
    site_pheno = stack(site_pheno, value_vars, ...
    'IndexVariableName', 'pheno', 'NewDataVariableName', 'doy');

    % Add site-level metadata to phenology data
    site_pheno{:, 'sitename'}        = data_meta.sitename(i);
    site_pheno{:, 'lat'}             = data_meta.lat(i);
    site_pheno{:, 'lon'}             = data_meta.lon(i);
    site_pheno{:, 'landcover_igbp'}  = data_meta.landcover_igbp(i);
    site_pheno{:, 'roi_id'}          = data_meta.roi_id(i);
    site_pheno{:, 'roi_veg_type'}    = data_meta.roi_veg_type(i);
    site_pheno{:, 'roi_n_maskfiles'} = data_meta.roi_n_maskfiles(i);

    cell_pheno{i} = site_pheno;
end

% Concatenate phenology data for all sites
data_pheno_raw_roi = vertcat(cell_pheno{:});


%% ------ STEP 6: Clean data (Outlier detection) ------ 
% Compute Z-scores for outlier detection
data_pheno_raw_roi = data_pheno_raw_roi(data_pheno_raw_roi.year >= year_start & ...
                                        data_pheno_raw_roi.year <= year_end, :);

data_pheno_raw = groupsummary(data_pheno_raw_roi, {'sitename', 'lat', 'lon', 'pheno', 'year'}, 'mean', 'doy'); 
data_pheno_raw.doy = data_pheno_raw.mean_doy;
data_pheno_raw(:, {'mean_doy', 'GroupCount'})= [];

stats = groupsummary(data_pheno_raw, {'sitename', 'lat', 'lon', 'pheno'}, {'mean', 'std'}, 'doy'); 
stats = stats(stats.GroupCount >= 3, :);

data_clean = innerjoin(data_pheno_raw, stats, 'Keys', {'sitename', 'lat', 'lon', 'pheno'});

data_clean.z_scores = (data_clean.doy - data_clean.mean_doy) ./ data_clean.std_doy;
data_clean          = data_clean(data_clean.z_scores <= z_score_thr, :);

data_phase          = table(phase_lvls', phase_lbls', 'VariableNames', ["pheno", "phase"]);
data_phase.pheno    = categorical(data_phase.pheno);
data_clean          = innerjoin(data_clean, data_phase, 'Keys', 'pheno');   
data_clean          = data_clean(:, {'sitename', 'lat', 'lon', 'phase', 'year', 'doy'});

writetable(data_clean, fullfile(dir_prep, 'data_ground_phenocam_clean.csv'));


%% ------ STEP 7: Spatial gridding ------
% Aggregate phenological data into grid cells based on latitude/longitude,
data_grid = data_clean;
[data_grid.row, data_grid.col] = get_rowcol_from_latlon(data_grid.lat, data_grid.lon, spatial_res);

% Apply circular mean for gridding (to handle the cyclic nature of DOY)
data_grid = groupsummary(data_grid, {'phase', 'row', 'col', 'year'}, ...
                         @(x) round(circ_mean_doy(x)), 'doy');
data_grid.Properties.VariableNames{'fun1_doy'} = 'doy';
data_grid.GroupCount = [];


%% ------ STEP 8: Plot for visual check ------
figure('Name', '[Pheno Cam] Greenup', 'Color', 'w');
histogram(data_grid.doy(data_grid.phase == "Greenup"));
xlabel('Day of Year (DOY)');
ylabel('Frequency');
title('[Pheno Cam] Greenup Distribution');

figure('Name', '[Pheno Cam] Dormancy', 'Color', 'w');
histogram(data_grid.doy(data_grid.phase == "Dormancy"));
xlabel('Day of Year (DOY)');
ylabel('Frequency');
title('[Pheno Cam] Dormancy Distribution')


%% ------ STEP 9: Save processed data ------
s_out = struct('data_grid', data_grid);
file_phenocam = sprintf('data_ground_phenocam_%d_%d.mat', min(data_grid.year), max(data_grid.year));
save(fullfile(dir_prep, file_phenocam), '-fromstruct', s_out, '-v7.3')