 % soildataQC.m
%
% Plots raw data and filtered soil temperature and moisture data 
% for a chosen SNOTEL site. 
% Performs statistics on these data and checks the filtering methods. 
% Useful for checking data loading and filtering procedures, determining 
% snowpack presence, and various filtering thresholds.
%
% Data for all years (available) is represented.
%
% This script was taken from (and should replace) site_ts_variability.m

close all; clear all;     % clear any figures and variables
fignum = 0;     % used to increment figure number for plots

% Add path to nanmean, etc.
addpath('~/data/code_resources/m_common/stat_tbx/');
% add path to test_loadsnotel.m
addpath('functiontests/');

% What filters to use?
filter1 = 'sigma';
threshold1 = 3;
filter2 = 'hampel';
threshold2 = 3;

% Ask user for site number
siteID = str2double(input('Which SNOTEL station?: ', 's'));

% Some data must be removed from the datasets provided by NRCS. These data
% fall into two categories:
%     1: Data that is clearly an instrument or human-caused error, perhaps
%        related to power supply, site maintenance, sensor failures, or 
%        data loss. These are usually identifiable by error code values,
%        level shifts, or other unusual patterns in the data.
%     2: Measurements that are outliers to the statistical distribution of 
%        the data. These are less frequent and reasons for them are hard to
%        pinpoint. Some could be power or instrument related, but they are
%        not always localized in time and are thus harder to identify and
%        remove.
% 
% The loadsnotel.m script is designed to catch the first category of error.
% There are a number of filtering functions that are designed to catch the
% other category. They are tested in the following blocks.

% First load hourly and daily data using both loadsnotel.m and a raw data
% loading function (test_loadsnotel.m).
hourlyData = loadsnotel('hourly', siteID);
dailyData = loadsnotel('daily', siteID);
hourlyRaw = test_loadsnotel('hourly', siteID);
dailyRaw = test_loadsnotel('daily', siteID);

% Parse out the date/times
decday_h = datenum(strcat(hourlyData{2}, hourlyData{3}), 'yyyy-mm-ddHH:MM');
years_h = floor(length(decday_h)/8760);
decday_d = datenum(dailyData{2}, 'yyyy-mm-dd');
decday_hRaw = datenum(strcat(hourlyRaw{2}, hourlyRaw{3}), 'yyyy-mm-ddHH:MM');
years_hRaw = floor(length(decday_hRaw)/8760);
decday_dRaw = datenum(dailyRaw{2}, 'yyyy-mm-dd');

% Snow water equivalent and airT
wteq = dailyData{4}; %nan*zeros(length(dailyData{4}), 23)];
airT = dailyData{9};
wteqRaw = dailyRaw{4}; %nan*zeros(length(dailyData{4}), 23)];
airTRaw = dailyRaw{9};
%wteq24hour = reshape(wteq', [], 1);

% Expand the daily swe data to hourly values by copying to a new array
% FIXME - Beware that time vectors don't match between raw and other 
% datasets, so a few numbers may be copied incorrectly (sub-daily in .csv
% files)
wteqHourly = zeros(length(decday_h), 1);
for i = 1:length(wteq)
%     index = decday_d == floor(decday_h(i));
%     wteqHourly(i) = wteq(index);
    wteqHourly(floor(decday_h)==decday_d(i)) = wteq(i);
end
wteqHourlyRaw = zeros(length(decday_hRaw), 1);
for i = 1:length(wteqRaw)
    index = decday_dRaw == floor(decday_hRaw(i));
    % Only use the first value because any second values are from
    % sub-daily rows.
%     values = wteqRaw(index);
%     wteqHourlyRaw(i) = values(1);
    wteqHourlyRaw(floor(decday_hRaw)==decday_dRaw(i)) = wteqRaw(i);
end

% Logical test for snowcover based on SWE values
swe_snowcover = wteqHourly > 0.5;
swe_snowcoverRaw = wteqHourlyRaw > 0.5;

% Assign some arrays for plotting
Ts5 = hourlyData{7}; % column 7 is at -2 in (5cm depth)
Ts20 = hourlyData{8}; % column 8 is at -8 in (20cm depth)
Ts50 = hourlyData{9}; % column 9 is at -20 in (50cm depth)
Ts5Raw = hourlyRaw{7}; % column 7 is at -2 in (5cm depth)
Ts20Raw = hourlyRaw{8}; % column 8 is at -8 in (20cm depth)
Ts50Raw = hourlyRaw{9}; % column 9 is at -20 in (50cm depth)
WC5 = hourlyData{4}; % column 4 is at -2 in (5cm depth)
WC20 = hourlyData{5}; % column 5 is at -8 in (20cm depth)
WC50 = hourlyData{6}; % column 6 is at -20 in (50cm depth)
WC5Raw = hourlyRaw{4}; % column 4 is at -2 in (5cm depth)
WC20Raw = hourlyRaw{5}; % column 5 is at -8 in (20cm depth)
WC50Raw = hourlyRaw{6}; % column 6 is at -20 in (50cm depth)

% TEST LOADDATA.M
% Plot raw data and data removed by loaddata.m (in red).
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', ['Site ' num2str(siteID) ' - Hourly Ts @ 3 depths']);
subplot(3, 1, 1)
plot(decday_hRaw, Ts5Raw, '.r', decday_h, Ts5, '.k');
title('Ts -5cm'); datetick('x', 'mmm-yy', 'keeplimits');
subplot(3, 1, 2)
plot(decday_hRaw, Ts20Raw, '.r', decday_h, Ts20, '.k');
title('Ts -20cm'); datetick('x', 'mmm-yy', 'keeplimits');
subplot(3, 1, 3)
plot(decday_hRaw, Ts50Raw, '.r', decday_h, Ts50, '.k');
title('Ts - 50cm'); datetick('x', 'mmm-yy', 'keeplimits')
legend('Raw data', 'After loadsnotel.m', 'Location', 'NorthWest');

fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', ['Site ' num2str(siteID) ' - Hourly VWC @ 3 depths']);
subplot(3, 1, 1)
plot(decday_hRaw, WC5Raw, '.r', decday_h, WC5, '.k');
title('VWC -5cm'); datetick('x', 'mmm-yy', 'keeplimits');
subplot(3, 1, 2)
plot(decday_hRaw, WC20Raw, '.r', decday_h, WC20, '.k');
title('VWC -20cm'); datetick('x', 'mmm-yy', 'keeplimits');
subplot(3, 1, 3);
plot(decday_hRaw, WC50Raw, '.r', decday_h, WC50, '.k');
title('VWC -50cm'); datetick('x', 'mmm-yy', 'keeplimits');
legend('Raw data', 'After loadsnotel.m', 'Location', 'NorthWest');

% TEST Ts FILTERING
% Generate filtered Ts data
 Ts5_filt1 = filtertempseries(Ts5, filter1, threshold1);
 Ts20_filt1 = filtertempseries(Ts20, filter1, threshold1);
 Ts50_filt1 = filtertempseries(Ts50, filter1, threshold1);
 Ts5_filt2 = filtertempseries(Ts5, filter2, threshold2);
 Ts20_filt2 = filtertempseries(Ts20, filter2, threshold2);
 Ts50_filt2 = filtertempseries(Ts50, filter2, threshold2);
% Ts_bothdiff = filterseries(Ts_meandiff, 'shift', 3);

% PLOT unfiltered and FILTER 1 timeseries
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', ['Site ' num2str(siteID) ' - Hourly Ts @ 3 depths']);
subplot(3, 1, 1)
plot(decday_h, Ts5, '.r', decday_h, Ts5_filt1, '.k');
title('Ts -5cm'); datetick('x', 'mmm-yy', 'keeplimits');
subplot(3, 1, 2)
plot(decday_h, Ts20, '.r', decday_h, Ts20_filt1, '.k');
title('Ts -20cm'); datetick('x', 'mmm-yy', 'keeplimits');
subplot(3, 1, 3)
plot(decday_h, Ts50, '.r', decday_h, Ts50_filt1, '.k');
title('Ts - 50cm'); datetick('x', 'mmm-yy', 'keeplimits');
legend('Raw data', ['Filtered using ' filter1], 'Location', 'NorthWest');

% PLOT histograms of unfiltered and FILTER 1 filtered data
fignum = fignum + 1;
h = figure(fignum);
set(h, 'Name', ['Site ' num2str(siteID) ' - ' ...
    ' Histogram of ' filter1 ' filter values']);
% Number of bins
xedges = linspace(-5, 5, 51);
%subplot(1,2,1);
n1 = histc(Ts5, xedges); % sort unfiltered values into bins
bar(xedges, n1, 'r'); % plot in barchart
hold on;
n2 = histc(Ts5_filt1, xedges);% sort mean diff values into bins
bar(xedges, n2, 'k');
ylabel({'number of','occurences'});
xlabel('Ts - 24hr mean Ts');

% PLOT unfiltered and FILTER 2 timeseries
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', ['Site ' num2str(siteID) ' - Hourly Ts @ 3 depths']);
subplot(3, 1, 1)
plot(decday_h, Ts5, '.r', decday_h, Ts5_filt2, '.k');
title('Ts -5cm'); datetick('x', 'mmm-yy', 'keeplimits');
subplot(3, 1, 2)
plot(decday_h, Ts20, '.r', decday_h, Ts20_filt2, '.k');
title('Ts -20cm'); datetick('x', 'mmm-yy', 'keeplimits');
subplot(3, 1, 3)
plot(decday_h, Ts50, '.r', decday_h, Ts50_filt2, '.k');
title('Ts - 50cm'); datetick('x', 'mmm-yy', 'keeplimits');
legend('Raw data', ['Filtered using ' filter2], 'Location', 'NorthWest');

% PLOT histograms of unfiltered and FILTER 2 data
fignum = fignum + 1;
h = figure(fignum);
set(h, 'Name', ['Site ' num2str(siteID) ' - ' ...
    ' Histogram of ' filter2 ' filter values']);
% Number of bins
xedges = linspace(-5, 5, 51);
%subplot(1,2,1);
n1 = histc(Ts5, xedges); % sort unfiltered values into bins
bar(xedges, n1, 'r'); % plot in barchart
hold on;
n2 = histc(Ts5_filt2, xedges);% sort mean diff values into bins
bar(xedges, n2, 'k');
ylabel({'number of','occurences'});
xlabel('Ts - shifted difference Ts');

% TEST VWC FILTERING
% Generate filtered VWC data
 WC5_filt1 = filtertempseries(WC5, filter1, threshold1);
 WC20_filt1 = filtertempseries(WC20, filter1, threshold1);
 WC50_filt1 = filtertempseries(WC50, filter1, threshold1);
 WC5_filt2 = filtertempseries(WC5, filter2, threshold2);
 WC20_filt2 = filtertempseries(WC20, filter2, threshold2);
 WC50_filt2 = filtertempseries(WC50, filter2, threshold2);
% Ts_bothdiff = filterseries(Ts_meandiff, 'shift', 3);

% PLOT unfiltered and filtered data
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', ['Site ' num2str(siteID) ' - Hourly VWC @ 3 depths']);
subplot(3, 1, 1)
plot(decday_h, WC5, '.r', decday_h, WC5_filt1, '.k');
title('VWC -5cm'); datetick('x', 'mmm-yy', 'keeplimits');
subplot(3, 1, 2)
plot(decday_h, WC20, '.r', decday_h, WC20_filt1, '.k');
title('VWC -20cm'); datetick('x', 'mmm-yy', 'keeplimits');
subplot(3, 1, 3)
plot(decday_h, WC50, '.r', decday_h, WC50_filt1, '.k');
title('VWC - 50cm'); datetick('x', 'mmm-yy', 'keeplimits');
legend('Raw data', ['Filtered using ' filter1], 'Location', 'NorthWest');

fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', ['Site ' num2str(siteID) ' - Hourly VWC @ 3 depths']);
subplot(3, 1, 1)
plot(decday_h, WC5, '.r', decday_h, WC5_filt2, '.k');
title('VWC -5cm'); datetick('x', 'mmm-yy', 'keeplimits');
subplot(3, 1, 2)
plot(decday_h, WC20, '.r', decday_h, WC20_filt2, '.k');
title('VWC -20cm'); datetick('x', 'mmm-yy', 'keeplimits');
subplot(3, 1, 3)
plot(decday_h, WC50, '.r', decday_h, WC50_filt2, '.k');
title('VWC - 50cm'); datetick('x', 'mmm-yy', 'keeplimits');
legend('Raw data', ['Filtered using ' filter2], 'Location', 'NorthWest');

% FIXME - this may no longer be necessary
% INTERPOLATION to fill gaps in Ts (helps with running mean and
% variance calculations below.
%
% Interpolates over NaNs (data gaps) in the input time series (may be
% complex), but ignores trailing and leading NaN.
%
% from FIXGAPS routine on Matlab Central file exchange,
% by R. Pawlowicz 6/Nov/99

Ts_filled = Ts5;

bad = isnan(Ts5);
good = find(~bad);

bad([1:(min(good)-1) (max(good)+1):end]) = 0;

Ts_filled(bad)=interp1(good, Ts5(good), find(bad), 'pchip');

Ts_filtered = Ts5_shiftdiff;

% PLOT original data over filled/filtered data to view differences
% fignum = fignum + 1;
% h = figure(fignum);
% set(h, 'Name', ['Site ' num2str(siteID) ' - Ts data filtering and interpolation']);
% subplot(211);
% plot(decday_h, Ts_filled, '.r');
% hold on
% plot(decday_h, Ts, '.k');
% title('Ts interpolation - Filled data points in red')
% 
% subplot(212);
% plot(decday_h, Ts_filled, '.r');
% hold on
% plot(decday_h, Ts_filtered, '.k');
% title('Ts filtering - Data from above in red, filtered to black series')

%If filled data looks good make Ts equal the filled data
Ts = Ts_filled;

% Create cell arrays of data and data labels
Ts_cell = {Ts, Ts_filtered};
Ts_label = {'UNFILTERED' 'SHIFT-DIFF'};

for i = 1:length(Ts_cell)
    
    Ts = Ts_cell{i};
    
    % COMPUTE a 24hour RUNNING mean & stddev
    % (For DISCRETE 24 hour intervals see version 1 of this program)
    window = 24; % 24 hour window
    % Calculate running mean
    runmean = filter(ones(window,1)/window, 1, Ts);
    % Calculate running variance - Should be something like
    % sum(xi-xbar).^2/(window-1);
    runvar = (filter(ones(window,1), 1, Ts.^2)-window*runmean.^2)/(window-1);
    
    % ERROR %
    % This variance calculation results in some negative values with this
    % dataset. This produces complex numbers in runstd below. For now filter
    % out these negative values with a logical test
    testvar = runvar < 0;
    runvar(testvar) = 0.009;
    
    % Shift back to  match data
    runmean = circshift(runmean, -0.5*window);
    runvar = circshift(runvar, -0.5*window);
    % Calculate running StdDev
    runstd = sqrt(runvar);
    % Another way to calculate moving standard deviation (doesn't work yet)
    %runstd2 = movingstd(Ts, 24);
    
    
    % FIND OUTLIERS in data (should differ b/t filtered and unfiltered Ts)
    % Difference between the running mean and each data point
    meandiff = Ts - runmean;
    % Difference between each datapoint and the datapoint one hour before it
    shiftdiff = Ts - circshift(Ts, 1);
    
    % PLOT Ts, SWE, AirT
    fignum = fignum + 1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - ' Ts_label{i} ...
        ' 5cm soil temperature, SWE, air temp']);
    subplot(3,1,1);
    plot(decday_h, Ts, ':r', decday_h, runmean, '-g', decday_h, runstd, 'b-');
    legend('5cm Ts', '24hr mean', '24hr StdDev');
    %
    % Ts and SWE
    subplot(3,1,2);
    [ax, ~, ~] = plotyy(decday_h, runmean, decday_d, wteq, 'plot');
    % Make more useful y axes
    set(ax(1), 'ylim', [-10 30]);
    %set(ax(1), 'ytick',  (0:0.2:5));
    set(ax(2), 'ylim', [0 40]);
    legend('24hr Mean Ts', 'SWE' );
    ylabel('24hr mean Ts.');
    %
    % StdDev of Ts and SWE on same plot
    subplot(3,1,3);
    plot(decday_d, airT, 'Color', [0.7 0.7 0.7]);
    hold on;
    plot(decday_h, runmean, 'r');
    ylim([-30 30]);
    legend('Air Temp', '24hr Mean Ts');
    ylabel('^oC');
    datetick('x', 'mmm-yy', 'keeplimits');
    
    % Comparison of 2009 and 2004 Mosby Mtn (AGU 2011 poster)
    % First set the dates to compare - a low snow and high snow year
    lowStart = datenum('Jul 1, 2009');
    lowEnd = lowStart + 300;
    highStart = datenum('Jul 1, 2005');
    highEnd = highStart + 300;
    testLow_d = decday_d > lowStart & decday_d < lowEnd;
    testHigh_d = decday_d > highStart & decday_d < highEnd;
    testLow_h = decday_h > lowStart & decday_h < lowEnd;
    testHigh_h = decday_h > highStart & decday_h < highEnd;
    
    fignum = fignum + 1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - ' Ts_label{i} ...
        ' 2 year comparison']);
    % Ts and SWE
    subplot(2,2,1);
    plot(decday_d(testHigh_d), wteq(testHigh_d), 'b', 'LineWidth', 1.5);
    xlim([highStart highEnd]);
    ylim([0 15]);
    ylabel('SWE (mm)');
    title('2005-6');
    %
    subplot(2,2,2);
    plot(decday_d(testLow_d), wteq(testLow_d), 'b');
    xlim([lowStart lowEnd]);
    ylim([0 15]);
    ylabel('SWE (mm)');
    title('2009-10');
    %
    subplot(2,2,3);
    plot(decday_d(testHigh_d), airT(testHigh_d), 'Color', [0.7,0.7,0.7],...
        'LineWidth', 1.5);
    hold on
    plot(decday_h(testHigh_h), runmean(testHigh_h), 'r');
    ylim([-25 20]);
    xlim([highStart highEnd]);
    % for some reason the axis changes randomly in the following lines so
    % be sure to set the xlimits (above)
    zeroline = line(get(gca, 'XLim'), [0, 0]);
    set(zeroline, 'Color', 'k', 'LineStyle', ':');
    legend('Air Temp', '24hr Mean Ts' );
    ylabel('Temp (^oC)');
    datetick('x','mmm', 'keeplimits');
    
    subplot(2,2,4);
    plot(decday_d(testLow_d), airT(testLow_d), 'Color', [0.7,0.7,0.7],...
        'LineWidth', 1.5);
    hold on
    plot(decday_h(testLow_h), runmean(testLow_h), 'r');
    ylim([-25 20]);
    xlim([lowStart lowEnd]);
    % for some reason the axis changes randomly in the following lines so
    % be sure to set the xlimits (above)
    zeroline = plot(get(gca, 'xlim'), [0, 0]);
    set(zeroline, 'Color', 'k', 'LineStyle', ':');
    %legend('Air Temp', '24hr Mean Ts' );
    ylabel('Temp (^oC)');
    datetick('x', 'mmm', 'keeplimits');
    
    
    % PLOT Ts, StdDev, outlier checks, and SWE (if available)
    fignum = fignum + 1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - ' Ts_label{i} ...
        ' 5cm soil temperature - rolling averages']);
    subplot(4,1,1);
    plot(decday_h, Ts, ':r', decday_h, runmean, '-g', decday_h, runstd, 'b-');
    legend('5cm Ts', '24hr mean', '24hr StdDev');
    %
    % Difference b/t Ts and mean (mean difference filter)
    subplot(4,1,2);
    plot(decday_h, meandiff, 'b-');
    ylabel('Degrees C');
    title('Running mean difference');
    %
    % Difference b/t Ts and adjacent Ts value (shift difference filter)
    subplot(4,1,3)
    plot(decday_h, shiftdiff, 'b-');
    ylabel('Degrees C');
    title('1hr shift difference');
    %
    % StdDev of Ts and SWE on same plot
    subplot(4,1,4);
    [ax, ~, ~] = plotyy(decday_h, runstd, decday_d, wteq, 'plot');
    % Make more useful y axes
    set(ax(1), 'ylim', [0 5]);
    set(ax(1), 'ytick',  (0:0.2:5));
    set(ax(2), 'ylim', [0 40]);
    legend('Ts 24hr StdDev', 'SWE' );
    ylabel('24hr rolling StdDev.');
    xlabel('Numeric day');

    % PLOT histograms of mean difference and shift difference filters
    fignum = fignum + 1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - ' Ts_label{i} ...
        ' Histograms of mean and shift filter values']);
    % Number of bins
    xedges = linspace(-5, 5, 51);
    subplot(1,2,1);
    n = histc(meandiff, xedges); % sort meandiff values into bins
    bar(xedges, n, 'k'); % plot in barchart
    ylabel({'number of','occurences'});
    xlabel('Ts - 24hr mean Ts');
    subplot(1,2,2);
    n = histc(shiftdiff, xedges);% sort shiftdiff values into bins
    bar(xedges, n, 'k');
    xlabel('Ts@t - Ts@t-1hr')


% By moving this endpoint down in file, other plots will run for UNFILTERED
% and FILTERED data
end;


% PLOT 24hr StdDev vs mean Ts (5cm)
fignum = fignum + 1;
h = figure(fignum);
set(h, 'Name', ['Site ' num2str(siteID) ' - ' Ts_label{i} ...
    ' 24hr SD vs Mean Ts']);
plot(runmean, runstd, 'g.')
title('Running');
xlabel('Mean Ts');
ylabel('24hr Running StdDev');


% SNOWCOVER THRESHOLDS and MEAN MONTHLY SOIL TEMPS
%
% Set the std deviation threshold for snowcover
sd_threshold = 0.1;
% Create a months array from the datenum array
tvec = datevec(decday_h);
months = tvec(:, 2);
% Use threshold to select snowcovered and snow free tests
% ERROR %
% the nans in runstd seem to be plotting so remove them with ~isnan test
snowtest = runstd < sd_threshold | isnan(runstd);
freetest = ~snowtest;

%PLOT timeseries of Ts during snowcovered and snowfree periods
fignum = fignum + 1;
h = figure(fignum);
set(h, 'Name', ['Site ' num2str(siteID) ' - ' Ts_label{i} ...
    ' Snowcover threshold tests']);
plot(decday_h(snowtest), Ts(snowtest), '-b', decday_h(freetest), ...
    Ts(freetest), '-r', decday_h, Ts, ':g', decday_h, runstd, 'm');
title(['Threshold = ' num2str(sd_threshold) ' degrees C, still produces snowfree points in winter']);
legend('Snowcover test', 'Snow-free test', 'All data');

% SUM up data for each month
sum_snowdays = zeros(12, 4);
sum_freedays = zeros(12, 4);
for j = 1:12
    monthtest = months == j;
    sum_snowdays(j,1) = j;
    sum_snowdays(j,2) = nanmean(Ts(snowtest & monthtest));
    sum_snowdays(j,3) = nanstd(Ts(snowtest & monthtest));
    sum_snowdays(j,4) = sum_snowdays(j,3)/sqrt(sum(snowtest & monthtest));
    sum_freedays(j,1) = j;
    sum_freedays(j,2) = nanmean(Ts(freetest & monthtest));
    sum_freedays(j,3) = nanstd(Ts(freetest & monthtest));
    sum_freedays(j,4) = sum_freedays(j,3)/sqrt(sum(freetest & monthtest));
end

% PLOT monthly snowcovered and snowfree temp + errors
fignum = fignum + 1;
h = figure(fignum);
set(h, 'Name', ['Site ' num2str(siteID) ' - ' Ts_label{i} ...
    ' Snow/snowfree mean Ts']);
errorbar(sum_snowdays(:,1), sum_snowdays(:,2), sum_snowdays(:,4));
hold on
errorbar(sum_freedays(:,1), sum_freedays(:,2), sum_freedays(:,4), 'r');
xlabel('Month');
ylabel('Celsius');
legend('Snow', 'Snow-free');
title(['Threshold = ' num2str(sd_threshold) ' degrees C']);

junk = 99;
