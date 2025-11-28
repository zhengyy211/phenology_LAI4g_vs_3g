function y_smoothed = w_savgol(y, varargin)
% -------------------------------------------------------------------------
% Function:    w_savgol
% Author:      Yaoyao Zheng
% Email:       zhengyaoyao@stu.pku.edu.cn
% Description: Apply a Savitzky-Golay filter to weighted and/or non-uniform spacing data.
%               This function performs Savitzky-Golay filtering with optional weights, non-uniform spacing,
%               and multiple iterations to update weights based on the filtered data.
%
% INPUTS:
%   y           - Original time series, [n, 1]
%   x           - (Optional) Index of y (for non-uniform spacing, default: 1:length(y))
%   weight      - (Optional) Weights for the data, default to ones([length(y), 1])
%   frame       - (Optional) Frame length (positive odd integer, default: 5)
%   polynom     - (Optional) Polynomial order (default: 2)
%   iter        - (Optional) Number of iterations (default: 1)
%
% OUTPUT:
%   y_smoothed  - Smoothed version of y after applying the Savitzky-Golay filter.
%
% Example usage:
%   y_smoothed = w_savgol(y);
%   y_smoothed = w_savgol(y, 'x', x); x is not the time of the time series, but the index
%   y_smoothed = w_savgol(y, 'weight', w);
%   y_smoothed = w_savgol(y, 'weight', w, 'x', x);
%   y_smoothed = w_savgol(y, 'weight', w, 'x', x, 'frame', 5, 'polynom', 2, 'iter', 3);
% 
% References:
%   1. Ranghetti, L. (2020). sen2rts: Build and Analyse Sentinel-2 Time Series.
%      https://rdrr.io/github/ranghetti/sen2rts/man/w_savgol.html
%   2. Chen, J. et al. (2004). A simple method for reconstructing a high-quality NDVI 
%      time-series data set based on the Savitzky–Golay filter. 
%      Remote Sensing of Environment. https://doi.org/10.1016/j.rse.2004.03.014.
%   3. Kong, D. et al. (2022). phenofit: An R package for extracting vegetation phenology from time series 
%      remote sensing. Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210X.13870
%   4. Nicholson, J. (2023). Savitzky-Golay Smoothing Filter. MATLAB Central File Exchange.
%      https://www.mathworks.com/matlabcentral/fileexchange/45420-savitzky-golay-smoothing-filter
% 
% Note:
%   This is a Matlab adaptation of R function at
%   https://rdrr.io/github/ranghetti/sen2rts/man/w_savgol.html,
%   with the addition of weight update function following Chen et al(2004).
% 
%   The R function is an adaptation of Python function at 
%   https://dsp.stackexchange.com/questions/1676/savitzky-golay-smoothing-filter-for-not-equally-spaced-data,
%   with the addition of weights following https://en.wikipedia.org/wiki/Savitzky-Golay_filter.
% 
%   For uniform spacing y without weights, this function will give the same output as sgolayfilt at
%   https://ww2.mathworks.cn/help/signal/ref/sgolayfilt.html
% -------------------------------------------------------------------------

% Initialize variables
p = inputParser;
addParameter(p, 'x', transpose(1:length(y)));
addParameter(p, 'weight', ones([length(y), 1]));
addParameter(p, 'frame', 5);
addParameter(p, 'polynom', 2);
addParameter(p, 'iter', 1);
parse(p, varargin{:});
x = p.Results.x;
weight = p.Results.weight;
frame = p.Results.frame;
polynom = p.Results.polynom;
iter = p.Results.iter;

% set weights of all NA values to 0
y = fillmissing(y, 'linear');
weight(isnan(y)) = 0;
weight(isnan(weight)) = 0;

if iter == 1
    y_smoothed = w_savgol_single(y, weight, x, frame, polynom);
else
    for ii = 1:1:iter
        y_smoothed = w_savgol_single(y, weight, x, frame, polynom);
        % https://doi.org/10.1016/j.rse.2004.03.014 
        % Section 2.2.3 - 2.2.5
        weight = wSG_weightUpdat(y, y_smoothed, weight);
        idx = y < y_smoothed;
        y(idx) = y_smoothed(idx);
    end
end

end


% A Matlab adaptation of R function at 
% https://rdrr.io/github/ranghetti/sen2rts/man/w_savgol.html
function y_smoothed = w_savgol_single(y, weight, x, frame, polynom)

half_frame = fix(frame / 2);
polynom = polynom + 1;

% Start smoothing
A = ones([frame, polynom]);
y_smoothed = zeros(length(y), 1);
for ii = (half_frame+1):1:(length(y)-half_frame)
    % Center a window of x values
    t = x(ii + (1:frame) - 1 - half_frame) - x(ii);

    % Diagonal matrix of weights
    w = weight((ii-half_frame):(ii+half_frame)); % vector of weights
    w = w * length(w) / sum(w, 'all'); % normalize
    W = diag(w); % diagonal matrix

    % Matrix A
    for j = 2:1:polynom
        A(:, j) = t .* A(:, j-1);
    end

    % Calculate c0 which is also the y value for y(ii)
    coeffs = (A.' * W * A) \ A.' * W;
    y_smoothed(ii) = coeffs(1, 1:frame) * y(ii + (1:frame) - 1 - half_frame);

    % If at the end or beginning, store all coefficients for the polynom
    if ii == half_frame+1
        first_coeffs = coeffs(1:polynom, 1:frame) * y(1:frame);
    elseif ii == (length(x) - half_frame - 1)
        last_coeffs  = coeffs(1:polynom, 1:frame) * y(length(y) - frame + (1:frame));
    end
end

% Interpolate the result at the left border
for ii = 1:1:half_frame
    x_i = 1;
    for jj = 1:1:polynom
        y_smoothed(ii) = y_smoothed(ii) + first_coeffs(jj) * x_i;
        x_i = x_i * (x(ii) - x(half_frame+1));
    end
end

% Interpolate the result at the right border
for ii = (length(y)-half_frame+1):1:length(y)
    x_i = 1;
    for jj = 1:1:polynom
      y_smoothed(ii) = y_smoothed(ii) + last_coeffs(jj) * x_i;
      x_i = x_i * (x(ii) - x(length(y) - half_frame));
    end
end

end


% The weight update function according to Chen et al (2004)
function wnew = wSG_weightUpdat(y, y_smoothed, weight)
  wnew = weight;
  re_abs = abs(y_smoothed - y);
  d_max = max(re_abs);
  I_pos = y_smoothed > y;
  wnew(I_pos) = (1 - re_abs(I_pos)./d_max);
end