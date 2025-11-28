function meta_table = read_phenocam_roi_meta(roi_file)
% -------------------------------------------------------------------------
% Function:    read_phenocam_roi_meta
% Author:      Yaoyao Zheng
% Email:       zhengyaoyao@stu.pku.edu.cn
% Description: This function reads metadata information from a Phenocam 
%              Region of Interest (ROI) file and extracts key details such 
%              as site name, ROI ID, vegetation type, and description. 
%              The extracted metadata is then returned as a table.
%
% INPUTS:
%   roi_file    - Path to the Phenocam ROI metadata file (string).
%
% OUTPUTS:
%   meta_table  - A table containing the metadata for the ROI, with the following fields:
%                 - roi_site_name:   Name of the site.
%                 - roi_id:          Identifier for the ROI.
%                 - roi_veg_type:    Type of vegetation for the ROI.
%                 - roi_description: Description of the ROI.
%                 - roi_n_maskfiles: Number of mask files associated with the ROI.
% -------------------------------------------------------------------------

% Open the ROI metadata file in read mode
fid = fopen(roi_file, 'r');
if fid == -1
    error('Cannot open file: %s', roi_file);
end

% Read the contents of the ROI metadata file into a table
roi_meta = readtable(roi_file, 'CommentStyle', '#');

% Get the number of mask files
n_maskfiles = height(roi_meta);

meta_struct = struct( ...
    'roi_site_name', "" , ...         % Name of the site
    'roi_id', NaN, ...                % ROI identifier (initialized to NaN)
    'roi_veg_type', "", ...           % Vegetation type (initialized to empty string)
    'roi_description', "", ...        % Description of the ROI (initialized to empty string)
    'roi_n_maskfiles', n_maskfiles);  % Number of mask files

% Process the file line by line to extract metadata
while ~feof(fid)
    line = strtrim(fgetl(fid));
    
    if startsWith(line, '# Site:')
        meta_struct.roi_site_name = string(strtrim(extractAfter(line, '# Site:')));
    elseif startsWith(line, '# ROI ID Number:')
        meta_struct.roi_id = str2double(strtrim(extractAfter(line, '# ROI ID Number:')));
    elseif startsWith(line, '# Veg Type:')
        meta_struct.roi_veg_type = string(strtrim(extractAfter(line, '# Veg Type:')));
    elseif startsWith(line, '# Description:')
        meta_struct.roi_description = string(strtrim(extractAfter(line, '# Description:')));
    end
end

fclose(fid);

meta_table = struct2table(meta_struct);

end