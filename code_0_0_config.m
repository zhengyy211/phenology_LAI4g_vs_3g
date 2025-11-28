%% ========================================================================
%  File:        code_0_0_config.m
%  Author:      Yaoyao Zheng
%  Email:       zhengyaoyao@stu.pku.edu.cn
%  MATLAB:      R2025b
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
cd(fullfile(cfg.path.root, 'scripts'));

% Add utility functions
addpath(fullfile(cfg.path.root, 'scripts', 'utils'));


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


%% ------ 3. Product Information ------
% LAI products and the time range for LAI
cfg.lai.products    = ["gimms_lai4g", "gimms_lai3g", "glass_lai", "yuan_lai"];
cfg.lai.year_start  = 1982;
cfg.lai.year_end    = 2015;

cfg.lai.gimms.doy_list = [  8,  23,  38,  52,  67,  82,  98, 113, 128, 143, 159, 174, ...
                          189, 204, 220, 235, 251, 266, 281, 296, 312, 327, 342, 357];
cfg.lai.glass_yuan.doy_list = 1:8:361;

% MODIS MCD12Q1 time range
cfg.modis.year_start = 2001;
cfg.modis.year_end   = 2015;

% Climate product
cfg.clim.year_start = cfg.lai.year_start - 1;
cfg.clim.year_end   = cfg.lai.year_end;
cfg.clim.doy_list   = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349];

% Define product-year combinations for processing
cfg.proc_list = {
    % {product, year_startart, year_end}
    "mcd12q2",           cfg.modis.year_start, cfg.modis.year_end;
    "yuan_lai",          cfg.modis.year_start, cfg.modis.year_end;
    cfg.lai.products(1), cfg.lai.year_start,   cfg.lai.year_end;
    cfg.lai.products(2), cfg.lai.year_start,   cfg.lai.year_end;  
    cfg.lai.products(1), cfg.lai.year_start,   cfg.modis.year_start - 1;
    cfg.lai.products(1), cfg.modis.year_start, cfg.modis.year_end;
    cfg.lai.products(2), cfg.lai.year_start,   cfg.modis.year_start - 1;
    cfg.lai.products(2), cfg.modis.year_start, cfg.modis.year_end;
    };


%% ------ 4. Spatial Configuration ------
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


%% ------ 5. Phenology Extraction Parameters ------
cfg.pheno.window_size_mth   = 3;  % Moving average window size (in months)
cfg.pheno.min_peak_dist_mth = 3;  % Minimum distance between peaks for detecting GS

% Growing season constraints
cfg.pheno.GS_min = 0.9;  % Minimum threshold for GS
cfg.pheno.GS_max = 1.1;  % Maximum threshold for GS

% Savitzky–Golay filter parameters
cfg.pheno.half_window = 2;  % Savitzky–Golay filter half-window
cfg.pheno.sg_framelen = 2 * cfg.pheno.half_window + 1; % Frame length for filtering

% Buffer size and outlier threshold
cfg.pheno.buffer_size = 31; % Days required beyond threshold to confirm SOS/EOS
cfg.pheno.z_score_thr = 3;  % Z-score threshold for outliers

% Shift in months for each season
cfg.pheno.season_shift_months = struct(...
    'spring',  3, ...   % spring -> shift +3 months
    'summer',  0, ...   % Summer -> no shift
    'autumn', -3, ...   % autumn -> shift -3 months
    'winter',  6 ...    % winter -> shift +6 months
);

% Define thresholds for different LAI phases (SOS, EOS, POS)
cfg.pheno.lai.trs_list = [0.15, 0.5];
trs_name = string(round(cfg.pheno.lai.trs_list*100));
cfg.pheno.lai.phases   = [append("SOS_", trs_name), "POS", append("EOS_", fliplr(trs_name))];
cfg.pheno.modis.phases = ["Greenup", "MidGreenup", "Peak", "MidGreendown", "Dormancy"];
clearvars trs_name

% Minimum number of in-situ observations required (e.g., from PEP725/USANPN)
cfg.pheno.min_obs = 3;

%% ------ 6. Sensitivity Analysis ------
% Variables for sensitivity analysis
cfg.sens.vars_clim = ["tmp", "prcp", "srad"];  % Temperature, Precipitation, Solar Radiation
cfg.sens.max_lag_mths = [1, 3:3:12]; % Maximum lag periods (in months) to consider for sensitivity analysis


%% ------ 7. Color Maps ------
% paletteer_d('RColorBrewer::PiYG')
% paletteer_d('RColorBrewer::PRGn')
% paletteer_d('RColorBrewer::BrBG')
% paletteer_d('RColorBrewer::RdBu')
cfg.cmp.PiYG = [ 39   100    25;
                 77   146    33;
                127   188    65;
                184   225   134;
                230   245   208;
                253   224   239;
                241   182   218;
                222   119   174;
                197    27   125;
                142     1    82]/255;
cfg.cmp.PRGn = [  0    68    27;
                 27   120    55;
                 90   174    97;
                166   219   160;
                217   240   211;
                231   212   232;
                194   165   207;
                153   112   171;
                118    42   131;
                 64     0    75]/255;
cfg.cmp.BrBG = [  0    60    48;
                  1   102    94;
                 53   151   143; 
                128   205   193; 
                199   234   229; 
                246   232   195;
                223   194   125;
                191   129    45;
                140    81    10;
                 84    48     5]/255;
cfg.cmp.RdBu = [  5    48    97
                 33   102   172;
                 67   147   195;
                146   197   222;
                209   229   240;
                253   219   199;
                244   165   130;
                214    96    77;
                178    24    43;
                103     0    31]/255;

%% ------ 8. Metadata ------
cfg.meta.author         = 'Yaoyao Zheng';
cfg.meta.mail           = 'zhengyaoyao@stu.pku.edu.cn';
cfg.meta.affiliation    = 'Peking University, Shenzhen Graduate School';
cfg.meta.matlab_version = version;

fprintf('[Info] Configuration initialized successfully.\n');