function meta_table = read_phenocam_site_meta(meta_file)
% -------------------------------------------------------------------------
% Function:    read_phenocam_site_meta
% Author:      Yaoyao Zheng
% Email:       zhengyaoyao@stu.pku.edu.cn
% Description: This function reads metadata information from a Phenocam site 
%              metadata file and extracts relevant details such as site name, 
%              geographic coordinates, elevation, start and end dates, 
%              vegetation type, and other attributes. The extracted metadata 
%              is then returned as a table.
%
% INPUTS:
%   meta_file   - Path to the Phenocam site metadata file (string).
%
% OUTPUTS:
%   meta_table  - A table containing the extracted metadata with the 
%                 following fields:
%                 - sitename: Name of the site.
%                 - long_name: Full site name or description.
%                 - lat: Latitude of the site.
%                 - lon: Longitude of the site.
%                 - elevation: Elevation of the site in meters.
%                 - date_start: Start date of the site data.
%                 - date_end: End date of the site data.
%                 - nimage: Number of images associated with the site.
%                 - primary_veg_type: Primary vegetation type.
%                 - dominant_species: Dominant species at the site.
%                 - landcover_igbp: Land cover classification based on IGBP.
%                 - flux_sitenames: Associated flux sites.
% -------------------------------------------------------------------------

% Open the metadata file for reading
fid = fopen(meta_file, 'r');
if fid == -1
    error('Cannot open file: %s', meta_file);
end

% Read the contents of the metadata file, line by line
meta_lines = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
meta_lines = meta_lines{1};

meta_struct = struct(...
    'sitename', "", ...          % Site name
    'long_name', "", ...         % Full site name
    'lat', NaN, ...              % Latitude
    'lon', NaN, ...              % Longitude 
    'elevation', NaN, ...        % Elevation 
    'date_start', NaT, ...       % Start date
    'date_end', NaT, ...         % End date 
    'nimage', NaN, ...           % Number of images
    'primary_veg_type', "", ...  % Primary vegetation type 
    'dominant_species', "", ...  % Dominant species
    'landcover_igbp', NaN, ...   % Land cover type 
    'flux_sitenames', "");       % Flux site names
 
% Loop through each line of the file to extract metadata
for i = 1:length(meta_lines)
    line = strtrim(meta_lines{i});

    if startsWith(line, 'sitename:')
        meta_struct.sitename = string(strtrim(extractAfter(line, 'sitename:')));
    elseif startsWith(line, 'long_name:')
        meta_struct.long_name = string(strtrim(extractAfter(line, 'long_name:')));
    elseif startsWith(line, 'lat:')
        meta_struct.lat = str2double(strtrim(extractAfter(line, 'lat:')));
    elseif startsWith(line, 'lon:')
        meta_struct.lon = str2double(strtrim(extractAfter(line, 'lon:')));
    elseif startsWith(line, 'elevation:')
        meta_struct.elevation = str2double(strtrim(extractAfter(line, 'elevation:')));
    elseif startsWith(line, 'date_start:')
        date_str = strtrim(extractAfter(line, 'date_start:'));
        meta_struct.date_start = datetime(date_str, 'InputFormat', 'yyyy-MM-dd', 'Format', 'yyyy-MM-dd');
    elseif startsWith(line, 'date_end:')
        date_str = strtrim(extractAfter(line, 'date_end:'));
        meta_struct.date_end = datetime(date_str, 'InputFormat', 'yyyy-MM-dd', 'Format', 'yyyy-MM-dd');
    elseif startsWith(line, 'nimage:')
        meta_struct.nimage = str2double(strtrim(extractAfter(line, 'nimage:')));
    elseif startsWith(line, 'primary_veg_type:')
        meta_struct.primary_veg_type = string(strtrim(extractAfter(line, 'primary_veg_type:')));
    elseif startsWith(line, 'dominant_species:')
        meta_struct.dominant_species = string(strtrim(extractAfter(line, 'dominant_species:')));
    elseif startsWith(line, 'landcover_igbp:')
        meta_struct.landcover_igbp = str2double(strtrim(extractAfter(line, 'landcover_igbp:')));
    elseif startsWith(line, 'flux_sitenames:')
        meta_struct.flux_sitenames = string(strtrim(extractAfter(line, 'flux_sitenames:')));
    end
end

meta_table = struct2table(meta_struct);

end