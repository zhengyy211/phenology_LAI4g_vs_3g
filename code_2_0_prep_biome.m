%% ========================================================================
%  File:    code_2_0_prep_biome.m
%  Author:  Yaoyao Zheng
%  Email:   zhengyaoyao@stu.pku.edu.cn
%  MATLAB:  R2026a
%  
%  Description: 
%  This script prepares MODIS MCD12Q1 data.
% =========================================================================

clear;

%% ------ STEP 1: Load configuration ------
run('code_0_0_config.m');

% Define directories for raw MODIS data and processed outputs
dir_raw    = fullfile(cfg.path.raw, 'MCD12Q1');
dir_prep   = cfg.path.prep;

% Define study period for MODIS
year_start = 2001;
year_end   = cfg.time.year_end;
year_list  = year_start:year_end;
n_years    = numel(year_list);

% Spatial grid parameters
n_rows     = cfg.space.n_rows;
n_cols     = cfg.space.n_cols;
R          = cfg.space.R;

types = ["Type1", "Type3"];

for type = types

    %% ------ STEP 2: Define Biome Types ------
    if strcmpi(type, 'Type1')
        % MCD12Q1 Land Cover Type 1: Annual International Geosphere-Biosphere Programme (IGBP) classification
        biome_types = table;
        biome_types.value = [    1;     2;     3;     4;     5;     6;     7;     8;     9;    10;    11;    12;    13;        14;    15;    16];
        biome_types.name  = ["ENF"; "EBF"; "DNF"; "DBF";  "MF"; "CSH"; "OSH"; "WSA"; "SAV"; "GRA"; "WET"; "CRO"; "UBL"; "CRO/NVM"; "SNO"; "NVG"];
        num_biomes = height(biome_types);
    else
        % MCD12Q1 Land Cover Type 3: Annual Leaf Area Index (LAI) classification
        biome_types = table;
        biome_types.value = [    1;     2;     3;     4;     5;     6;     7;     8;     9;    10];
        biome_types.name  = ["GRA"; "SHR"; "CRO"; "SAV"; "EBF"; "DBF"; "ENF"; "DNF"; "NVG"; "UBL"];
        num_biomes = height(biome_types);
    end

    %% ------ STEP 3: Load the Land Cover Data ------
    cell_biome = cell(n_years, 1);

    for idx_yr = 1:n_years

        current_year = year_list(idx_yr);
        fprintf('[Info] Year %d...\n', current_year);

        %% ------ STEP 2.1: Read East/West tiles ------
        file_name_W = sprintf('MCD12Q1_LC_%s_%04d_W.tif', type, current_year);
        file_name_E = sprintf('MCD12Q1_LC_%s_%04d_E.tif', type, current_year);

        data_W = imread(fullfile(dir_raw, file_name_W));
        data_E = imread(fullfile(dir_raw, file_name_E));

        temp_biome = single([data_W, data_E]);
        temp_biome(temp_biome == 255) = NaN;
        temp_biome(~ismember(temp_biome, biome_types.value)) = NaN;

        cell_biome{idx_yr} = imresize_block_mode(temp_biome, n_rows, n_cols);

    end

    data_biome = cat(3, cell_biome{:});
    data_biome = mode(data_biome, 3);
    [n_rows, n_cols] = size(data_biome, [1 2]);

    %% ------ STEP 4: Save Processed Data (MAT) ------
    s_out = struct;
    s_out.biome_types = biome_types;
    s_out.num_biomes  = num_biomes;
    s_out.data_biome  = reshape(data_biome, n_rows * n_cols, 1);

    path_out = fullfile(dir_prep, sprintf('data_biome_%s.mat', type));
    save(path_out, '-fromstruct', s_out, '-v7.3')

    geotiffwrite(fullfile(dir_prep, sprintf('map_biome_%s.tif', type)), data_biome, R);

end