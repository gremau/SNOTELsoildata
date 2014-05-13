% soildataQC.m
%
% Plots raw data and filtered soil temperature and moisture data 
% for a chosen SNOTEL site. Also checks some of the numbers determined by
% summarize_wateryear.m (snowpack onset, melt, etc).
% 
% Useful for checking data loading and filtering procedures, snowpack 
% presence/abscense, and various filtering thresholds.
%
% Data for all years (available) is represented.
%
% This script was taken from (and replaces) site_ts_variability.m

close all; clear all;     % clear any figures and variables
fignum = 0;     % used to increment figure number for plots

% Add paths to nanmean, hline/vline etc.
addpath('~/data/code_resources/m_common/stat_tbx/');
addpath('~/data/code_resources/m_common/hline_vline/');
% add path to loadsnotel_raw.m
addpath('testing/');

% Get a list of hourly files
hrFilelist = csvread('../processed_data/filelist_hourly.txt');
hrFilelist = unique(hrFilelist(:,1));

% Ask user for site number
siteID = str2double(input('Which SNOTEL station?: ', 's'));

% Ask user whether to run loadsnotel() in 'exclude mode', removing bad data
% in the exclude files.
exclude = input('Exclude bad data with loadsnotel.m?  (y/n) : ', 's');

% What filter parameters to use?
filter1 = 'sigma';
threshold1 = 3;
% Mean filter sucks, stopped using it
%filter2 = 'mean';
%threshold2 = 3;

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
% other category. All are tested in the following blocks.

% First load hourly and daily data using both loadsnotel.m and a raw data
% loading function (loadsnotel_raw.m).
if strcmpi(exclude, 'y') && any(hrFilelist==siteID)
    hourlyData = loadsnotel(siteID, 'hourly', 'exclude');
    hourlyRaw = loadsnotel_raw(siteID, 'hourly');
    dailyData = loadsnotel(siteID, 'daily', 'exclude');
elseif strcmpi(exclude, 'y') && ~any(hrFilelist==siteID)
    disp('Hourly file not found');
    hourlyData = cell(1, 10);
    hourlyRaw = cell(1, 10);
    dailyData = loadsnotel(siteID, 'daily', 'exclude');
elseif ~strcmpi(exclude, 'y') && any(hrFilelist==siteID)
    hourlyData = loadsnotel(siteID, 'hourly');
    hourlyRaw = loadsnotel_raw(siteID, 'hourly');
    dailyData = loadsnotel(siteID, 'daily');
elseif ~strcmpi(exclude, 'y') && ~any(hrFilelist==siteID)
    disp('Hourly file not found');
    hourlyData = cell(1, 10);
    hourlyRaw = cell(1, 10);
    dailyData = loadsnotel(siteID, 'daily');
end
dailyRaw = loadsnotel_raw(siteID, 'daily');

% Load the climate data from summarize_wateryear.m
processeddatapath = '../processed_data/';
climData = csvread([processeddatapath 'wyear_climatesummary.txt'], 1,0);

% Parse out the date/times
decday_h = datenum(strcat(hourlyData{2}, hourlyData{3}), 'yyyy-mm-ddHH:MM');
wyears_h = hourlyData{10};
%years_h = floor(length(decday_h)/8760);
decday_d = datenum(dailyData{2}, 'yyyy-mm-dd');
wyears_d = dailyData{21};
decday_hRaw = datenum(strcat(hourlyRaw{2}, hourlyRaw{3}), 'yyyy-mm-ddHH:MM');
%years_hRaw = floor(length(decday_hRaw)/8760);
decday_dRaw = datenum(dailyRaw{2}, 'yyyy-mm-dd');

% Assign some arrays for plotting
wteq = dailyData{4}; % SWE;
precip = dailyData{5}; % Precip
snowd = dailyData{10}; % Snow depth
airTobs = dailyData{6}; % Air temps
airTmax = dailyData{7};
airTmin = dailyData{8};
airTavg = dailyData{9};  
wteqRaw = dailyRaw{4}; % SWE;
precipRaw = dailyRaw{5}; % Precip
snowdRaw = dailyRaw{10}; % Snow depth
airTobsRaw = dailyRaw{6}; % Air temps
airTmaxRaw = dailyRaw{7};
airTminRaw = dailyRaw{8};
airTavgRaw = dailyRaw{9};

ts5 = hourlyData{7}; % column 7 is at -2 in (5cm depth)
ts20 = hourlyData{8}; % column 8 is at -8 in (20cm depth)
ts50 = hourlyData{9}; % column 9 is at -20 in (50cm depth)
ts5Raw = hourlyRaw{7}; % column 7 is at -2 in (5cm depth)
ts20Raw = hourlyRaw{8}; % column 8 is at -8 in (20cm depth)
ts50Raw = hourlyRaw{9}; % column 9 is at -20 in (50cm depth)
vwc5 = hourlyData{4}; % column 4 is at -2 in (5cm depth)
vwc20 = hourlyData{5}; % column 5 is at -8 in (20cm depth)
vwc50 = hourlyData{6}; % column 6 is at -20 in (50cm depth)
vwc5Raw = hourlyRaw{4}; % column 4 is at -2 in (5cm depth)
vwc20Raw = hourlyRaw{5}; % column 5 is at -8 in (20cm depth)
vwc50Raw = hourlyRaw{6}; % column 6 is at -20 in (50cm depth)

% Parse out climate data for the site
climData = climData(climData(:, 1)==siteID, :);
minWyDatenums = datenum(climData(:,2)-1, 10, 1);
maxSWE = climData(:,3);
maxSWEDatenums = minWyDatenums(:) + climData(:,4);
onsetDatenums = minWyDatenums(:) + climData(:,6);
meltDatenums = minWyDatenums(:) + climData(:,7);

% Expand the daily swe data to hourly boolean values using swe_snowcover.m
wteqHourlyBool = swe_snowcover(siteID, decday_h);

% ******* TEST LOADDATA.M *******
% Plot raw data and data removed by loadsnotel.m (in red).
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', ['Site ' num2str(siteID) ' - Daily PRECIP data']);
subplot(3, 1, 1)
plot(decday_dRaw, wteqRaw, '.r', decday_d, wteq, '.k');
title('SWE (in)'); datetick('x', 'mmm-yy');
subplot(3, 1, 2)
plot(decday_dRaw, precipRaw, '.r', decday_d, precip, '.k');
title('Precip (in)'); datetick('x', 'mmm-yy');
subplot(3, 1, 3)
plot(decday_dRaw, snowdRaw, '.r', decday_d, snowd, '.k');
title('SnowDepth (in)'); datetick('x', 'mmm-yy')
legend('Raw data', 'After loadsnotel.m', 'Location', 'NorthWest');

fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', ['Site ' num2str(siteID) ' - Daily AIR T data']);
subplot(4, 1, 1)
plot(decday_dRaw, airTobsRaw, '.r', decday_d, airTobs, '.k');
title('Observed AirT'); datetick('x', 'mmm-yy');
subplot(4, 1, 2)
plot(decday_dRaw, airTmaxRaw, '.r', decday_d, airTmax, '.k');
title('Max AirT'); datetick('x', 'mmm-yy');
subplot(4, 1, 3)
plot(decday_dRaw, airTminRaw, '.r', decday_d, airTmin, '.k');
title('Min AirT'); datetick('x', 'mmm-yy');
subplot(4, 1, 4)
plot(decday_dRaw, airTavgRaw, '.r', decday_d, airTavg, '.k');
title('Avg AirT'); datetick('x', 'mmm-yy')
legend('Raw data', 'After loadsnotel.m', 'Location', 'NorthWest');

if any(hrFilelist==siteID)
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - Hourly Ts @ 3 depths']);
    subplot(3, 1, 1)
    plot(decday_hRaw, ts5Raw, '.r', decday_h, ts5, '.k');
    title('Ts -5cm'); datetick('x', 'mmm-yy');
    subplot(3, 1, 2)
    plot(decday_hRaw, ts20Raw, '.r', decday_h, ts20, '.k');
    title('Ts -20cm'); datetick('x', 'mmm-yy');
    subplot(3, 1, 3)
    plot(decday_hRaw, ts50Raw, '.r', decday_h, ts50, '.k');
    title('Ts - 50cm'); datetick('x', 'mmm-yy')
    legend('Raw data', 'After loadsnotel.m', 'Location', 'NorthWest');
    
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - Hourly VWC @ 3 depths']);
    subplot(3, 1, 1)
    plot(decday_hRaw, vwc5Raw, '.r', decday_h, vwc5, '.k');
    title('VWC -5cm'); datetick('x', 'mmm-yy');
    subplot(3, 1, 2)
    plot(decday_hRaw, vwc20Raw, '.r', decday_h, vwc20, '.k');
    title('VWC -20cm'); datetick('x', 'mmm-yy');
    subplot(3, 1, 3);
    plot(decday_hRaw, vwc50Raw, '.r', decday_h, vwc50, '.k');
    title('VWC -50cm'); datetick('x', 'mmm-yy');
    legend('Raw data', 'After loadsnotel.m', 'Location', 'NorthWest');
    
    % ******* TEST Ts FILTERING *******
    % Generate filtered Ts data
    ts5_F1 = filterseries(ts5, filter1, 25, threshold1);
    ts20_F1 = filterseries(ts20, filter1, 25, threshold1);
    ts50_F1 = filterseries(ts50, filter1, 25, threshold1);
    
    % PLOT unfiltered and FILTER 1 timeseries
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - Hourly Ts @ 3 depths']);
    subplot(3, 1, 1)
    plot(decday_h, ts5, '.r', decday_h, ts5_F1, '.k');
    title('Ts -5cm'); datetick('x', 'mmm-yy', 'keeplimits');
    subplot(3, 1, 2)
    plot(decday_h, ts20, '.r', decday_h, ts20_F1, '.k');
    title('Ts -20cm'); datetick('x', 'mmm-yy', 'keeplimits');
    subplot(3, 1, 3)
    plot(decday_h, ts50, '.r', decday_h, ts50_F1, '.k');
    title('Ts - 50cm'); datetick('x', 'mmm-yy', 'keeplimits');
    legend('Raw data', ['Filtered using ' filter1], 'Location', 'NorthWest');
    
    %******* TEST VWC FILTERING *******
    % Generate filtered VWC data
    vwc5_F1 = filterseries(vwc5, filter1, 25, threshold1);
    vwc20_F1 = filterseries(vwc20, filter1, 25, threshold1);
    vwc50_F1 = filterseries(vwc50, filter1, 25, threshold1);
    
    % PLOT unfiltered and FILTER 1 timeseries
    fignum = fignum+1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - Hourly VWC @ 3 depths']);
    subplot(3, 1, 1)
    plot(decday_h, vwc5, '.r', decday_h, vwc5_F1, '.k');
    title('VWC -5cm'); datetick('x', 'mmm-yy', 'keeplimits');
    subplot(3, 1, 2)
    plot(decday_h, vwc20, '.r', decday_h, vwc20_F1, '.k');
    title('VWC -20cm'); datetick('x', 'mmm-yy', 'keeplimits');
    subplot(3, 1, 3)
    plot(decday_h, vwc50, '.r', decday_h, vwc50_F1, '.k');
    title('VWC - 50cm'); datetick('x', 'mmm-yy', 'keeplimits');
    legend('Raw data', ['Filtered using ' filter1], 'Location', 'NorthWest');
    
    %****HISTOGRAMS of unfiltered and filtered Tsoil and VWC data******
    fignum = fignum + 1;
    h = figure(fignum);
    set(h, 'Name', ['Site ' num2str(siteID) ' - ' ...
        ' Histogram of filtered values']);
    % Number of bins
    xedges_ts = linspace(-20, 47.5, 51);
    xedges_vwc = linspace(-1, 47.5, 51);
    % Sort the raw and filtered series into bins
    uts = histc(ts5, xedges_ts);    % Unfiltered data
    uvwc = histc(vwc5, xedges_vwc);    % Unfiltered data
    fts = histc(ts5_F1, xedges_ts); % Filter 1
    fvwc = histc(vwc5_F1, xedges_vwc); % Filter 2
    % Now plot the histograms
    subplot(1,2,1);
    bar(xedges_ts, uts, 'r'); % Unfiltered values in red
    hold on;
    bar(xedges_ts, fts, 'k'); % Filtered values in black
    title(['Filtered with ' filter1]);
    ylabel({'Frequency'}); xlabel('Tsoil');
    subplot(1,2,2);
    bar(xedges_vwc, uvwc, 'r'); % Unfiltered values in red
    hold on;
    bar(xedges_vwc, fvwc, 'k'); % Filtered values in black
    title(['Filtered with ' filter1]);
    ylabel('Frequency'); xlabel('VWC');
end

% ******* TEST PRECIP and AIR T FILTERING *******
wteq_F1 = filterseries(wteq, filter1, 11, threshold1);
precip_F1 = filterseries(precip, filter1, 11, threshold1);
snowd_F1 = filterseries(snowd, filter1, 11, threshold1);
airTobs_F1 = filterseries(airTobs, filter1, 11, threshold1);
airTmax_F1 = filterseries(airTmax, filter1, 11, threshold1);
airTmin_F1 = filterseries(airTmin, filter1, 11, threshold1);
airTavg_F1 = filterseries(airTavg, filter1, 11, threshold1);

% PLOT unfiltered and FILTER 1 precip timeseries
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', ['Site ' num2str(siteID) ' - PRECIP Filter 1']);
subplot(3, 1, 1)
plot(decday_d, wteq , '.r', decday_d, wteq_F1, '.k');
title('Ts -5cm'); datetick('x', 'mmm-yy');
subplot(3, 1, 2)
plot(decday_d, precip, '.r', decday_d, precip_F1, '.k');
title('Ts -20cm'); datetick('x', 'mmm-yy');
subplot(3, 1, 3)
plot(decday_d, snowd, '.r', decday_d, snowd_F1, '.k');
title('Snow depth'); datetick('x', 'mmm-yy');
legend('From loadsnotel.m', ['Filtered using ' filter1], 'Location', 'NorthWest');

% PLOT unfiltered and FILTER 1 airT timeseries
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', ['Site ' num2str(siteID) ' - AIR T Filter 1']);
subplot(4, 1, 1)
plot(decday_d, airTobs, '.r', decday_d, airTobs_F1, '.k');
title('Observed AirT'); datetick('x', 'mmm-yy');
subplot(4, 1, 2)
plot(decday_d, airTmax, '.r', decday_d, airTmax_F1, '.k');
title('Max AirT'); datetick('x', 'mmm-yy');
subplot(4, 1, 3)
plot(decday_d, airTmin, '.r', decday_d, airTmin_F1, '.k');
title('Min AirT'); datetick('x', 'mmm-yy');
subplot(4, 1, 4)
plot(decday_d, airTavg, '.r', decday_d, airTavg_F1, '.k');
title('Avg AirT'); datetick('x', 'mmm-yy')
legend('from loadsnotel.m', ['Filtered using ' filter1], 'Location', 'NorthWest');


% ******* TEST HOURLY SWE and ONSET/MELT CALCULATIONS ********

% PLOT Ts, SWE, AirT
fignum = fignum + 1;
h = figure(fignum);
set(h, 'Name', ['Site ' num2str(siteID) ' - ' ...
    ' Ts, SWE, and onset/melt calculations']);
subplot(3,1,1);
plot(decday_d, wteq, '-b', decday_h, wteqHourlyBool*15, '--g');
datetick('x', 'mmm-yy');
legend('SWE', 'Hourly snowcover boolean (swe_snowcover.m)');
% SWE, peak SWE, onset, melt (as calculated in summarize_wateryear.m)
subplot(3,1,2);
plot(decday_d, wteq, '-b');
hold on;
plot(maxSWEDatenums, maxSWE, '*r');
vline(onsetDatenums, ':b');
vline(meltDatenums, ':k');
legend('SWE', 'Peak SWE', 'Onset day', 'Snowmelt day'); 
datetick('x', 'mmm-yy');
% Ts, with onset and melt timing
subplot(3,1,3);
plot(decday_h, ts5_F1, '-r');
hold on;
vline(onsetDatenums, ':b');
vline(meltDatenums, ':k');
legend('Ts - 5cm'); 
datetick('x', 'mmm-yy');

% Some example sites

if siteID==762
    
    fignum = fignum+1;
    
    fig1 = figure('position',[100 0 800 450],'paperpositionmode',...
        'auto', 'color','none','InvertHardcopy','off');
    set(fig1, 'Name', ['Data exclusion and filtering for site 762'])
    set(fig1, 'DefaultAxesFontSize',11, 'DefaultTextFontSize', 11);
    
    subplot(2, 1, 1)
    dlerrors = ts50Raw < -90;
    incwy = decday_hRaw > datenum('Oct 1 2008') & ...
        decday_hRaw < datenum('Sept 30 2009')
    h1 = plot(decday_hRaw, ts50Raw, '.r');
    hold on
    h2 = plot(decday_h, ts50, '.k');
    h3 = plot(decday_hRaw(dlerrors), ts50Raw(dlerrors), '.b');
    h4 = plot(decday_hRaw(incwy), ts50Raw(incwy), '.g');
    ylim([-105, 55]);
    text(0.95, 0.9, 'a','Units', 'normalized');
    datetick('x', 'mmm-yy', 'keeplimits');
    legend([h1, h3, h4], 'Data removed after visual inspection',...
        'Datalogger error', 'Incomplete water year',...
        'Location', 'southeast');
    set(gca, 'position', [1 .9 1 1.1] .* get(gca, 'position'));
    
    subplot(2, 1, 2)
    h1 = plot(decday_h, ts50, '.r');
    hold on
    h2 = plot(decday_h, ts50_F1, '.k');
    xlim([datenum('Sept 30, 2006'), datenum('Oct 1, 2011')]);
    text(0.95, 0.9, 'b','Units', 'normalized');
    datetick('x', 'mmm-yy', 'keeplimits');
    legend([h1], 'Filtered data (> 3 S.D.)', 'Location', 'south');
    text(-.07, .93,'T_{soil} (^oC)','Units', 'normalized', 'rotation', 90);
    set(gca, 'position', [1 1 1 1.1] .* get(gca, 'position'));
    
    figpath = '../figures/';
    print(fig1,'-depsc2','-painters',[figpath 'tsoil_exclude.eps']);
    figpath = '../../manuscript_1/figs/';
    print(fig1,'-depsc2','-painters',[figpath 'tsoil_exclude.eps'])
    
elseif siteID==550
    
    fignum = fignum+1;
    
    fig2 = figure('position',[100 0 800 450],'paperpositionmode',...
        'auto', 'color','none','InvertHardcopy','off');
    set(fig2, 'Name', ['Data exclusion and filtering for site 550'])
    set(fig2, 'DefaultAxesFontSize',11, 'DefaultTextFontSize', 11);
    
    subplot(2, 1, 1)
    dlerrors = vwc50Raw < -90;
    %incwy = decday_hRaw > datenum('Oct 1 2008') & ...
    %    decday_hRaw < datenum('Sept 30 2009')
    h1 = plot(decday_hRaw, vwc50Raw, '.r');
    hold on
    h2 = plot(decday_h, vwc50, '.k');
    h3 = plot(decday_hRaw(dlerrors), vwc50Raw(dlerrors), '.b');
    %h4 = plot(decday_hRaw(incwy), vwc50Raw(incwy), '.g');
    ylim([-105, 55]);
    text(0.95, 0.9, 'a','Units', 'normalized');
    datetick('x', 'mmm-yy', 'keeplimits');
    legend([h1, h3], 'Data removed after visual inspection',...
        'Datalogger error',...
        'Location', 'southeast');
    set(gca, 'position', [1 .9 1 1.1] .* get(gca, 'position'));
    
    subplot(2, 1, 2)
    h1 = plot(decday_h, vwc50, '.r');
    hold on
    h2 = plot(decday_h, vwc50_F1, '.k');
    xlim([datenum('Sept 30, 2006'), datenum('Oct 1, 2011')]);
    text(0.95, 0.9, 'b','Units', 'normalized');
    datetick('x', 'mmm-yy', 'keeplimits');
    legend([h1], 'Filtered data (> 3 S.D.)', 'Location', 'south');
    text(-.07, .93,'\theta (%)','Units', 'normalized', 'rotation', 90);
    set(gca, 'position', [1 1 1 1.1] .* get(gca, 'position'));
    
    figpath = '../figures/';
    print(fig2,'-depsc2','-painters',[figpath 'theta_exclude.eps']);
    figpath = '../../manuscript_1/figs/';
    print(fig2,'-depsc2','-painters',[figpath 'theta_exclude.eps'])
    
elseif siteID==342
    
    fignum = fignum+1;
    
    fig3 = figure('position',[100 0 800 250],'paperpositionmode',...
        'auto', 'color','none','InvertHardcopy','off');
    set(fig3, 'Name', ['Calibration problem at site 342'])
    set(fig3, 'DefaultAxesFontSize',11, 'DefaultTextFontSize', 11);
    
    dlerrors = vwc50Raw < -90;
    %incwy = decday_hRaw > datenum('Oct 1 2008') & ...
    %    decday_hRaw < datenum('Sept 30 2009')
    h1 = plot(decday_hRaw, vwc50Raw, '.r');
    hold on
    h2 = plot(decday_h, vwc50, '.k');
    %h3 = plot(decday_hRaw(dlerrors), vwc50Raw(dlerrors), '.b');
    %h4 = plot(decday_hRaw(incwy), vwc50Raw(incwy), '.g');
    ylim([-5, 55]);
    datetick('x', 'mmm-yy', 'keeplimits');
    legend(h1, 'Data removed after visual inspection',...
        'Location', 'northwest');
    ylabel('\theta (%)');
    %set(gca, 'position', [1 .9 1 1.1] .* get(gca, 'position'));
    
    figpath = '../figures/';
    print(fig3,'-depsc2','-painters',[figpath 'calproblem.eps']);
    figpath = '../../manuscript_1/figs/';
    print(fig3,'-depsc2','-painters',[figpath 'calproblem.eps'])
    
end



junk = 99;
