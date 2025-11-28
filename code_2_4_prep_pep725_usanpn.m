%% ========================================================================
%  File:        code_2_4_prep_pep725_usanpn.m
%  Author:      Yaoyao Zheng
%  Email:       zhengyaoyao@stu.pku.edu.cn
%  MATLAB:      R2025b
%  Description: This script reads, filters, cleans, and spatially grids the
%               phenological observation data from the PEP725 and USANPN networks.
%  Reference:   Zohner et al., Science, 2023. DOI:10.1126/science.adf5098
% =========================================================================

clear;

%% ------ STEP 1: Load configuration ------
run('code_0_0_config.m');

% dir_pep725   = fullfile(cfg.path.raw, 'pep725');
% dir_usanpn   = fullfile(cfg.path.raw, 'usanpn');
dir_pep725   = 'I:/PEP725';
dir_usanpn   = 'I:/USANPN';
dir_prep     = cfg.path.prep;

year_start   = cfg.lai.year_start;
year_end     = cfg.lai.year_end;
year_range   = year_start:year_end;

spatial_res  = cfg.space.spatial_res;
lat_min      = cfg.space.lat_limit(1);

min_obs      = cfg.pheno.min_obs;

s_raw = struct;

% Define spring and autumn phenophases for each dataset
s_raw.pep725.phase_spring = 11; % First leaves unfolded (BBCH code)
s_raw.pep725.phase_autumn = 95; % 50% of leaves fallen (BBCH code)

s_raw.usanpn.event_spring = [
    "Breaking leaf buds", ...
    "Breaking leaf buds (lilac/honeysuckle)", ...
    "Breaking needle buds (conifers)", ...
    "Breaking needle buds (deciduous)", ...
    "Initial growth (forbs)", ...
    "Initial growth (grasses/sedges)", ...
    "Emerging needles (pines)"];
s_raw.usanpn.event_autumn = [
    ">=50% of leaves fallen (deciduous)", ...
    ">=50% of needles fallen (larch)", ...
    "Falling leaves", ...
    "Falling needles"];

%% ------ STEP 2: Read and filter PEP725 data ------
%  Objective:
%       (1) Load all available PEP725 CSV files
%       (2) Retain only necessary variables (mapped via 'vars_map_pep725')
%       (3) Filter observations by target years and BBCH phenophase codes

vars_map_pep725 = {
    's_id',        's_id'; ...
    'lon',         'lon'; ...
    'lat',         'lat'; ...
    'genus',       'genus'; ...
    'species',     'species'; ...
    'gss_id',      'plant_id'; ...
    'phase_id',    'phase_id'; ...
    'qc_flag',     'qc_flag'; ...
    'cult_season', 'cult_season'; ...
    'year',        'year'; ...
    'day',         'doy'};

BBCH_codes = [s_raw.pep725.phase_spring, s_raw.pep725.phase_autumn];

file_list = dir(fullfile(dir_pep725, 'pep725_*.csv'));
n_file    = length(file_list);

cell_pep725 = cell(n_file, 1);
for idx = 1:n_file

    file_path = fullfile(file_list(idx).folder, file_list(idx).name);

    opts = detectImportOptions(file_path, 'Delimiter', ';');
    opts.SelectedVariableNames = vars_map_pep725(:,1);

    T = readtable(file_path, opts);
    for i = 1:size(vars_map_pep725, 1)
        T.Properties.VariableNames{i} = vars_map_pep725{i, 2};
    end

    % Filter based on year and BBCH codes
    T = T(ismember(T.year, year_range), :);
    T = T(ismember(T.phase_id, BBCH_codes), :);
    T = T(T.lat >= lat_min, :);

    cell_pep725{idx} = T;
end

s_raw.pep725_raw = vertcat(cell_pep725{:});


%% ------ STEP 3: Read and filter USANPN data ------
%  Objective:
%       (1) Load the USANPN "Leaves" and "Needles" datasets
%       (2) Harmonize variable names
%       (3) Identify target spring and autumn phenophases based on description
%       (4) Filter by temporal range and relevant phenophase IDs
vars_map_usanpn = { ...
    'Site_ID',                       's_id'; ...
    'Longitude',                     'lon'; ...
    'Latitude',                      'lat'; ...
    'Genus',                         'genus'; ...
    'Species',                       'species'; ...
    'Individual_ID',                 'plant_id'; ...
    'Phenophase_ID',                 'phase_id'; ...
    'Phenophase_Description',        'phase_desc'; ...
    'Phenophase_Category',           'phase_cat'; ...
    'First_Yes_Year',                'year'; ...
    'First_Yes_DOY',                 'doy'; ...
    'NumDays_Since_Prior_No',        'numDays_since_prior_no'; ...
    'NumDays_Until_Next_No',         'numDays_Until_next_no'; ...
    'Multiple_Observers',            'multiple_observers'; ...
    'Observed_Status_Conflict_Flag', 'observed_status_conflict_flag'};

data_leaves  = readtable(fullfile(dir_usanpn, 'Leaves',  'individual_phenometrics_data.csv'));
data_needles = readtable(fullfile(dir_usanpn, 'Needles', 'individual_phenometrics_data.csv'));
data_usanpn_raw = vertcat(data_leaves, data_needles);

% Subset and rename columns
data_usanpn_raw = data_usanpn_raw(:, vars_map_usanpn(:, 1));
for i = 1:size(vars_map_usanpn,1)
    data_usanpn_raw.Properties.VariableNames{i} = vars_map_usanpn{i,2};
end

% Convert categorical/text columns to string type for robust matching
data_usanpn_raw.phase_desc = string(data_usanpn_raw.phase_desc);
data_usanpn_raw.observed_status_conflict_flag = string(data_usanpn_raw.observed_status_conflict_flag);

% Derive phase IDs corresponding to spring and autumn descriptions
s_raw.usanpn.phase_spring = unique(data_usanpn_raw{ismember(data_usanpn_raw.phase_desc, s_raw.usanpn.event_spring), 'phase_id'});
s_raw.usanpn.phase_autumn = unique(data_usanpn_raw{ismember(data_usanpn_raw.phase_desc, s_raw.usanpn.event_autumn), 'phase_id'});
phase_ids = [s_raw.usanpn.phase_spring; s_raw.usanpn.phase_autumn];

% Temporal and phenological filtering
data_usanpn_raw = data_usanpn_raw(ismember(data_usanpn_raw.year, year_range), :);
data_usanpn_raw = data_usanpn_raw(ismember(data_usanpn_raw.phase_id, phase_ids), :);
data_usanpn_raw = data_usanpn_raw(data_usanpn_raw.lat >= lat_min, :);

s_raw.usanpn_raw = data_usanpn_raw;


%% ------ STEP 4: Data cleaning and spatial gridding ------
%  Objective:
%       For each dataset (PEP725 and USANPN):
%       (1) Apply dataset-specific quality control (QC) filters
%       (2) Clean phenological records
%       (3) Aggregate observations spatially by percentile binning
%       (4) Generate DOY histograms for visual QC

for prod = ["pep725", "usanpn"]

    fprintf('[Info] Cleaning %s data...\n', prod);

    phase_spring = s_raw.(prod).phase_spring;
    phase_autumn = s_raw.(prod).phase_autumn;
    data_raw     = s_raw.(sprintf('%s_raw', prod));

    % ------ STEP 4.1: Quality control filtering ------
    if prod == "pep725"
        % Retain only QC-verified (flag = 0) records
        data_qc = data_raw(data_raw.qc_flag == 0, :);
    else
        % Remove inconsistent or multi-observer conflicts
        data_qc = data_raw(~ismember(data_raw.observed_status_conflict_flag, ["MultiObserver-StatusConflict", "OneObserver-StatusConflict"]), :);
        data_qc = data_qc(data_qc.numDays_since_prior_no > 0 & data_qc.numDays_since_prior_no <= 14, :);
    end

    % ------ STEP 4.2: Phenological data cleaning ------
    group_vars = {'s_id', 'lon', 'lat', 'plant_id', 'phase_id'};
    phase_col  = 'phase_id';
    data_clean = clean_pheno_obs(data_qc, group_vars, phase_col, phase_spring, phase_autumn, min_obs);

    data_clean.phase = strings(height(data_clean),1);
    data_clean.phase(ismember(data_clean.phase_id, phase_spring)) = "Spring";
    data_clean.phase(ismember(data_clean.phase_id, phase_autumn)) = "Autumn";

    file_clean = sprintf('data_ground_%s.csv', prod);
    writetable(data_clean, fullfile(dir_prep, file_clean));

    % ------ STEP 4.3: Spatial gridding ------
    % Aggregate phenological data into grid cells based on latitude/longitude,
    % using percentile statistics (30th percentile threshold by default)
    data_grid = data_clean;
    [data_grid.row, data_grid.col] = get_rowcol_from_latlon(data_grid.lat, data_grid.lon, spatial_res);
  
    data_grid = groupsummary(data_grid, {'phase', 'row', 'col', 'year'}, @(x) round(prctile(x, 30)), 'doy');
    data_grid.Properties.VariableNames{'fun1_doy'} = 'doy';
    data_grid.GroupCount = [];


    % ------ STEP 4.4: Plot for visual check ------
    figure('Name', sprintf('[%s] Spring', prod), 'Color', 'w');
    histogram(data_grid.doy(data_grid.phase == "Spring"));
    xlabel('Day of Year (DOY)');
    ylabel('Frequency');
    title(sprintf('[%s] Spring Distribution', prod));

    figure('Name', sprintf('[%s] Autumn', prod), 'Color', 'w');
    histogram(data_grid.doy(data_grid.phase == "Autumn"));
    xlabel('Day of Year (DOY)');
    ylabel('Frequency');
    title(sprintf('[%s] Autumn Distribution', prod));

    %% ------ STEP 5: Save processed phenological data ------
    %  The final structure 's_pheno_obs' contains all intermediate and final
    %  datasets, including:
    %      - raw    (original filtered data)
    %      - clean  (after QC and phenophase separation)
    %      - grid   (spatially aggregated)

    s_pheno = struct('data_grid', data_grid);
    file_pheno = sprintf('data_ground_%s_%d_%d.mat', prod, min(data_grid.year), max(data_grid.year));
    save(fullfile(dir_prep, file_pheno), '-fromstruct', s_pheno, '-v7.3')
end

