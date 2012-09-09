function array = filtertempseries(series, type, threshold)
% filtertemprseries.m
%
% Filters soil temperature data series using a difference from the mean, or 
% difference from neighbor, then fills in generated nans with interpolation
% routine.
%
% arg 1 = input data series
% arg 2 = 'mean', 'median', 'shift', 'sigma', or 'hampel' filter type
% arg 3 = threshold difference above which datapoint is set to nan
%
% Note that if NaN's are present in the input series, most filters will
% them by the window size. Therefore, these calculate statistics based on
% interpolated data, generated using interpseries(series).
%
% 4/27/2011 Greg Maurer     

addpath('~/data/code_resources/m_common/slidefun/');
addpath('~/data/code_resources/m_common/movingstd/');
addpath('~/data/code_resources/m_common/hampel/');

% Moving window size - best to make this an ODD number.
window = 25; 

% MEAN - Filter by difference from the mean
% WARNING !! This filter propagates NaN's if interpolation is left out.
if strcmp(type, 'mean')
    % First interpolate over the missing data in the input series
    series_filled = interpseries(series);
    % Then calculate a running mean with the window size
    runningMean = filter(ones(window,1)/window, 1, series_filled);
    % Resulting mean is shifted forward in phase by window/2, shift it back
    runningMean = circshift(runningMean, -floor(window/2));
    % Find difference from the mean
    diff = series - runningMean;
    % Change datapoints more than the threshold value away from the 
    % mean to nan
    filteredSeries = series;
    testDiff = abs(diff) > threshold;
    filteredSeries(testDiff) = nan;

% MEDIAN - Filter by difference from the median
% Use slidefun.m from MATLAB FEx to calculate median
% WARNING !! This filter propagates NaN's if interpolation is left out.
elseif strcmp(type, 'median')
    % First interpolate over the missing data in the input series
    series_filled = interpseries(series);
    % Then calculate a running median based on the window size.
    runningMedian = slidefun(@median, window, series_filled);
    % Find difference from the mean
    diff = series - runningMedian;
    % Change datapoints more than the threshold value away from the 
    % mean to nan
    filteredSeries = series;
    testDiff = abs(diff) > threshold;
    filteredSeries(testDiff) = nan;
    
% SHIFT - Filter by difference from nearest neighbor;
% Warning - multiplies the NaN's in the original data
elseif strcmp(type, 'shift')
    %Calculate difference from nearest neigboring datapoint
    diff1 = series - circshift(series, 1);
    diff2 = series - circshift(series, -1);
    % Change datapoints more than the threshold value away from the 
    % neighbor to nan
    filteredSeries = series;
    testDiff = abs(diff1) > threshold | abs(diff2) > threshold;
    filteredSeries(testDiff) = nan;

% SIGMA - Filter that removes data using standard deviation of the data.
% Use movingstd.m from MATLAB FEx to calculate StdDev of series
% FIXME - Seems to remove a ton of winter Ts data for some reason - 
% investigate. Also - consider using a running median rather than mean 
% here, see Ron Pearson's ideas about this.
elseif strcmp(type, 'sigma')
    % First interpolate over the missing data in the input series
    series_filled = interpseries(series);
    % Then calculate a running std based on the window size.
    runningStd = movingstd(series_filled, window);
    % Also calculate a running mean and a sigmas vector
    runningMean = filter(ones(window,1)/window, 1, series_filled);
    sigmas = threshold * runningStd;
    % Resulting mean is shifted forward in phase by window/2, shift it back
    runningMean = circshift(runningMean, -floor(window/2));
    % Set a hi and low value (# of sigmas from mean) to filter outliers.
    hi = runningMean + sigmas;
    lo = runningMean - sigmas;
    % Change datapoints beyond the hi/lo thresholds to nan
    filteredSeries = series;
    testDiff = filteredSeries > hi | filteredSeries < lo;
    filteredSeries(testDiff) = nan;

% HAMPEL - Filter using a Hampel filter
% Use hampel.m from MATLAB FEx to filter outliers. This is a complicated
% algorithm, so look at the documentation before changing parameters.
elseif strcmp(type, 'hampel')
    x = 1:length(series)';
    % Hampel filter - parameters are default here, including a window size
    % that is based on the size of the input array, and a threshold value
    % of 3 (for removing outliers).
    [filteredSeries,remove,~,~,~,~,~] = hampel(x, series);
    filteredSeries(remove) = nan;
    
else
    error('Invalid filter type (mean, median, shift, sigma, or hampel)')
end

array = filteredSeries;