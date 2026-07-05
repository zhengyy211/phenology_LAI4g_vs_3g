%% ========================================================================
%  File:    code_2_3_prep_pep725.m
%  Author:  Yaoyao Zheng
%  Email:   zhengyaoyao@stu.pku.edu.cn
%  MATLAB:  R2026a
%
%  Description:
%  This script reads, filters, and cleans the phenological observation data
%  from the PEP725 networks.
%
%  Reference:   Zohner et al., Science, 2023. DOI:10.1126/science.adf5098
% =========================================================================

clear;

%% ------ STEP 1: Load configuration ------
run('code_0_0_config.m');

dir_raw      = fullfile(cfg.path.raw, 'pep725');
dir_prep     = cfg.path.prep;

year_start   = cfg.time.year_start;
year_end     = cfg.time.year_end;
year_range   = year_start:year_end;
year_suffix  = sprintf('%d_%d', year_start, year_end);

spatial_res  = cfg.space.spatial_res;
lat_min      = cfg.space.lat_limit(1);

min_obs      = cfg.pheno.min_obs;

% Define spring and autumn phenophases for each dataset
phase_spring = [10, 11, 182]; % BBCH code
phase_autumn = 95;            % 50% of leaves fallen (BBCH code)
phase_all    = [phase_spring, phase_autumn];

%% ------ STEP 2: Read and filter data ------
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
    'cult_season', 'cult_season'; ...
    'year',        'year'; ...
    'day',         'doy'};

file_list = dir(fullfile(dir_raw, 'pep725_*.csv'));
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

    % Filter based on lat, year and phase
    T = T(T.lat  >= lat_min & T.year >= year_start, :);
    T = T(ismember(T.phase_id, phase_all), :);

    cell_pep725{idx} = T;
end

pep725_raw         = vertcat(cell_pep725{:});
pep725_raw         = pep725_raw(ismember(pep725_raw.year, year_range), :);
pep725_raw.genus   = string(pep725_raw.genus);
pep725_raw.species = string(pep725_raw.species);
pep725_raw.species = regexprep(pep725_raw.species, '\s+', ' ');
pep725_raw         = pep725_raw(pep725_raw.cult_season ~= 2, :);

file_raw = sprintf('data_ground_pep725_raw_%s.csv', year_suffix);
writetable(pep725_raw, fullfile(dir_prep, file_raw));

%% ------ STEP 3: Data cleaning ------
group_vars   = {'s_id', 'lon', 'lat', 'genus', 'species', 'plant_id', 'phase_id'};
phase_col    = 'phase_id';
pep725_clean = clean_pheno_obs(pep725_raw, group_vars, phase_col, phase_spring, phase_autumn, min_obs);

pep725_clean.phase = strings(height(pep725_clean),1);
pep725_clean.phase(ismember(pep725_clean.phase_id, phase_spring)) = "Spring";
pep725_clean.phase(ismember(pep725_clean.phase_id, phase_autumn)) = "Autumn";

file_clean  = sprintf('data_ground_pep725_clean_%s.csv', year_suffix);
writetable(pep725_clean, fullfile(dir_prep, file_clean));

%% ------ STEP 4: Spatial gridding ------
data_grid = pep725_clean;
[data_grid.row, data_grid.col] = get_rowcol_from_latlon(data_grid.lat, data_grid.lon, spatial_res);

data_grid = groupsummary(data_grid, {'phase', 'row', 'col', 'year'}, @(x) round(prctile(x, 30)), 'doy');
data_grid.Properties.VariableNames{'fun1_doy'} = 'doy';
data_grid.GroupCount = [];

figure('Name', 'Spring', 'Color', 'w');
histogram(data_grid.doy(data_grid.phase == "Spring"));
xlabel('Day of Year (DOY)');
ylabel('Frequency');
title('Spring');

figure('Name', 'Autumn', 'Color', 'w');
histogram(data_grid.doy(data_grid.phase == "Autumn"));
xlabel('Day of Year (DOY)');
ylabel('Frequency');
title('Autumn');

%% ------ STEP 5: Save processed phenological data ------
file_grid = sprintf('data_ground_pep725_%s.mat', year_suffix);
save(fullfile(dir_prep, file_grid), 'data_grid', '-v7.3')