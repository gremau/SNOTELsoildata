function array = filtertempseries(series, type, threshold)
% filtertemprseries.m
%
% Filters soil temperature data series using a difference from the mean, or 
% difference from neighbor, then fills in generated nans with interpolation
% routine.
%
% arg 1 = input data series
% arg 2 = 'mean' or 'shift' difference filter type
% arg 3 = threshold difference above which datapoint is set to nan
%
% v1: 4/27/2011 Greg Maurer     

addpath('~/data/code_resources/m_common/');

% Filter by difference from the mean
if strcmp(type, 'mean')
    % First calculate a running mean
    window = 24;
    runningMean = filter(ones(window,1)/window, 1, series);
    % Shift back to  match data
    runningMean = circshift(runningMean, -0.5*window);
    % Find difference from the mean
    diff = series - runningMean;
    % Change datapoints more than the threshold value away from the 
    % mean to nan
    filteredSeries = series;
    testDiff = abs(diff) > threshold;
    filteredSeries(testDiff) = nan;
    
% Filter by difference from nearest neighbor;    
elseif strcmp(type, 'shift')
    %Calculate difference from nearest neigboring datapoint
    diff1 = series - circshift(series, 1);
    diff2 = series - circshift(series, -1);
    % Change datapoints more than the threshold value away from the 
    % neighbor to nan
    filteredSeries = series;
    testDiff = abs(diff1) > threshold | abs(diff2) > threshold;
    filteredSeries(testDiff) = nan;
else
    error('Not a valid filter type (mean or shift)')
end


% INTERPOLATION to fill gaps in tne filtered series (helps with running 
%mean and variance calculations used in other routines.
%
% Interpolates over NaNs (data gaps) in the input time series (may be
% complex), but ignores trailing and leading NaN.
%
% from FIXGAPS routine on Matlab Central file exchange,
% by R. Pawlowicz 6/Nov/99

filteredFilledSeries = filteredSeries;

bad = isnan(filteredSeries);
good = find(~bad);

bad([1:(min(good)-1) (max(good)+1):end]) = 0;

filteredFilledSeries(bad)=interp1(good, filteredSeries(good), find(bad), 'pchip');

% PLOT original data over interpolated data to view differences
% h = figure;
% set(h, 'Name', ['Filterseries.m data interpolation']);
% plot(series, '.r');
% hold on
% plot(runningMean, '-g');
% plot(filteredFilledSeries, '.b');
% plot(filteredSeries, '.k');
% title('Filtered points in red, interpolated data in blue')

% If filled data looks good return filtered and interpolated data
array = filteredFilledSeries;
