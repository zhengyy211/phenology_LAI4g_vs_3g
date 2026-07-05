%% ========================================================================
%  File:    code_0_0_config.m
%  Author:  Yaoyao Zheng
%  Email:   zhengyaoyao@stu.pku.edu.cn
%  MATLAB:  R2026a
%  
%  Description: 
%  Configuration file for reproducing the experiments described in the paper:
%
%  "Spatiotemporally Consistent Long-Term Remote Sensing Time Series 
%   Improve Phenology Monitoring in the Northern Hemisphere"
%
%  This script sets up the required directory paths, spatial-temporal settings,
%  and phenological extraction parameters needed for reproducing the experiments.
% =========================================================================


%% ------ 1. Environment Initialization ------
cfg = struct;

% Project root directory
cfg.path.root = 'J:/pheno_3g4g';
assert(exist(cfg.path.root,'dir') ~= 0, 'Root directory not found')
cd(cfg.path.root);

% Add utility functions
addpath(fullfile(cfg.path.root, 'scripts'))
addpath(fullfile(cfg.path.root, 'scripts', 'utils'))


%% ------ 2. Directory Structure ------
% Define the paths for different data directories
cfg.path.data   = fullfile(cfg.path.root, 'data');      % Root directory for all data
cfg.path.raw    = fullfile(cfg.path.data, '01_raw');    % Raw data (downloaded from various sources)
cfg.path.prep   = fullfile(cfg.path.data, '02_prep');   % Preprocessed data
cfg.path.pheno  = fullfile(cfg.path.data, '03_pheno');  % Extracted phenology
cfg.path.sens   = fullfile(cfg.path.data, '04_sens');   % Sensitivity analysis
cfg.path.valid  = fullfile(cfg.path.data, '05_valid');  % Validation of phenology
cfg.path.robust = fullfile(cfg.path.data, '06_robust'); % Data and results for robustness analysis
cfg.path.figure = fullfile(cfg.path.root, 'figure');    % Directory for figures

% Create directories if not exist
path_list = struct2cell(cfg.path);
for idx = 1:numel(path_list)
    if ~exist(path_list{idx}, 'dir')
        mkdir(path_list{idx});
    end
end
clearvars path_list idx


%% ------ 3. Spatial Temporal Configuration ------
cfg.time.year_start = 1982;
cfg.time.year_mid   = 2000;
cfg.time.year_end   = 2015;
cfg.time.year_end_l = 2024;

% Set latitude and longitude limits for the Northern Hemisphere
cfg.space.lat_limit   = [  30,  90];  % Northern Hemisphere (≥ 30°N)
cfg.space.lon_limit   = [-180, 180];
cfg.space.spatial_res = 1/12;         % Spatial resolution in degrees

% Calculate grid size based on latitude and longitude limits
cfg.space.n_rows = (cfg.space.lat_limit(2) - cfg.space.lat_limit(1)) / cfg.space.spatial_res;
cfg.space.n_cols = (cfg.space.lon_limit(2) - cfg.space.lon_limit(1)) / cfg.space.spatial_res;

% Geographic reference
cfg.space.R = georefcells(cfg.space.lat_limit, cfg.space.lon_limit, ...
    [cfg.space.n_rows, cfg.space.n_cols], ...
    'ColumnsStartFrom', 'north', 'RowsStartFrom', 'west');

% Find row indices corresponding to latitude range
lat_vals = ((90 - cfg.space.spatial_res/2):-cfg.space.spatial_res:(-90 + cfg.space.spatial_res/2));
cfg.space.row_st = find(lat_vals <= cfg.space.lat_limit(2), 1, 'first');
cfg.space.row_ed = find(lat_vals >= cfg.space.lat_limit(1), 1,  'last');
clearvars lat_vals


%% ------ 4. Product ------
%% ------ 4.1 VI ------
cfg.vi = struct;
cfg.vi.prods_gimms_lai = ["gimms_lai4g",  "gimms_lai3g"];
cfg.vi.prods_ohters    = ["glass_lai", "yuan_lai", "gimms_ndvi4g", "gimms_ndvi3g", "mc_ndvi"];
cfg.vi.products        = [cfg.vi.prods_gimms_lai, cfg.vi.prods_ohters];

for prod = cfg.vi.products

    % prod = cfg.vi.products(1)
    if contains(prod, 'modis') || contains(prod, 'yuan')
        year_start = max(cfg.time.year_start, 2001);
    else
        year_start = max(cfg.time.year_start, 1982);
    end

    switch lower(prod)
        case {'gimms_lai4g', 'gimms_lai3g', 'yuan_lai', 'gimms_ndvi4g'}
            year_end = cfg.time.year_end;
        case 'glass_lai'
            year_end = min(cfg.time.year_end_l, 2018);
        case 'gimms_ndvi3g'
            year_end = min(cfg.time.year_end_l, 2022);
        case 'mc_ndvi'
            year_end = min(cfg.time.year_end_l, 2025);
    end

    if contains(prod, 'glass') || contains(prod, 'yuan')
        doy_list = 1:8:366;
    elseif strcmpi(prod, 'mc_ndvi')
        doy_list = 1:16:366;
    else
        doy_list = [  7,  23,  38,  52,  66,  82,  97, 112, 127, 143, 158, 173, ...
                    188, 204, 219, 235, 250, 265, 280, 296, 311, 326, 341, 357];
    end

    if contains(prod, 'gimms')
        spatial_res = 1/12; 
    elseif contains(prod, 'glass') || strcmpi(prod, 'mc_ndvi')
        spatial_res = 0.05;
    elseif contains(prod, 'yuan')
        spatial_res = 0.1;
    end

    cfg.vi.(prod) = struct( ...
        'name',         replace(upper(prod), '_', ' '), ...
        'year_start',   year_start, ...
        'year_end',     year_end, ...
        'doy_list',     doy_list, ...
        'spatial_res',  spatial_res);

    clearvars year_start year_end doy_list spatial_res
end

% MCD12Q2
cfg.vi.mcd12q2 = struct( ...
    'name',          'MCD12Q2', ...
    'year_start',    2001, ...
    'year_end',      cfg.time.year_end);

% MOD15A2H
cfg.vi.modis_lai = struct( ...
    'name',          'MOD15A2H', ...
    'year_start',    2001, ...
    'year_end',      cfg.time.year_end, ...
    'spatial_res',   1/12, ... % 500 m -> 1/12° (mean)
    'doy_list',      1:8:361);

products = [cfg.vi.products, 'mcd12q2'];
n_prods  = numel(products);
proc_list = cell(n_prods, 3);
for ii = 1:n_prods
    prod = products(ii);
    proc_list(ii, :) = {prod, cfg.vi.(prod).year_start, cfg.vi.(prod).year_end};
end
proc_list = cell2table(proc_list, 'VariableNames', {'prod', 'year_start', 'year_end'});

proc_list_p1 = proc_list(ismember(proc_list.prod, cfg.vi.prods_gimms_lai), :);
proc_list_p1{:, 'year_end'}   = cfg.time.year_mid;
proc_list_p2 = proc_list(ismember(proc_list.prod, cfg.vi.prods_gimms_lai), :);
proc_list_p2{:, 'year_start'} = cfg.time.year_mid + 1;
proc_list = [proc_list; proc_list_p1; proc_list_p2];
cfg.vi.prod_list = proc_list;
clearvars products n_prods proc_list proc_list_p1 proc_list_p2 prod ii

%% ------ 4.2 Climate ------
cfg.clim.year_start = cfg.time.year_start - 1;
cfg.clim.year_end   = min(cfg.time.year_end_l, 2024);
cfg.clim.doy_list   = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349];


%% ------ 5. Phenology Extraction Parameters ------
cfg.pheno.window_size_mth   = 3;  % Moving average window size (in months)
cfg.pheno.min_peak_dist_mth = 3;  % Minimum distance between peaks for detecting GS

% Growing season constraints
cfg.pheno.GS_min = 0.9;  % Minimum threshold for GS
cfg.pheno.GS_max = 1.1;  % Maximum threshold for GS

% Buffer size and outlier threshold
cfg.pheno.buffer_size = 31; % Days required beyond threshold to confirm SOS/EOS
cfg.pheno.z_score_thr = 3;  % Z-score threshold for outliers

% Define thresholds for different phases (SOS, EOS, POS)
cfg.pheno.trs_list = [0.15, 0.5];
trs_name = string(round(cfg.pheno.trs_list*100));
cfg.pheno.phase_lvls = [append("SOS_", trs_name), "POS", append("EOS_", fliplr(trs_name))];
cfg.pheno.phase_lbls = ["Greenup", "MidGreenup", "Peak", "MidGreendown", "Dormancy"];
clearvars trs_name

% Minimum number of in-situ observations required
cfg.pheno.min_obs = 10;


%% ------ 6. Sensitivity Analysis ------
% Variables for sensitivity analysis
cfg.sens.vars_clim    = ["temp", "prcp", "srad"];  % Temperature, Precipitation, Solar Radiation
cfg.sens.max_lag_mths = [1, 3:3:12]; % Maximum lag periods (in months) to consider for sensitivity analysis


%% ------ 7. Metadata ------
cfg.meta.author         = 'Yaoyao Zheng';
cfg.meta.mail           = 'zhengyaoyao@stu.pku.edu.cn';
cfg.meta.affiliation    = 'Peking University, Shenzhen Graduate School';
cfg.meta.matlab_version = '26.1.0.3276743 (R2026a) Update 3';

fprintf('[Info] Configuration initialized successfully.\n');