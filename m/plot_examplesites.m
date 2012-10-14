% plot_examplesites.m
close all;
clear all;

% -----------------------------------------------------------------------
% Compare a site with a low and high early season snowpack
% Load Mosby Mtn hourly and daily data
siteID = 643; % 643=Mosby mtn, 
hourlyData = loadsnotel(siteID, 'hourly');
dailyData = loadsnotel(siteID, 'daily', 'exclude');

% Parse out the date/times
decday_h = datenum(strcat(hourlyData{2}, hourlyData{3}), 'yyyy-mm-ddHH:MM');
decday_d = datenum(dailyData{2}, 'yyyy-mm-dd');

% Parse out some variables
wteq = dailyData{4} * 25.4;
airT = dailyData{9};
ts = filterseries(hourlyData{7}, 'sigma', 25, 3); %7 = 5cm, 8 = 20cm, etc

% Comparison of 2009 and 2004 Mosby Mtn (AGU 2011 poster)
% First set the dates to compare - a low snow and high snow year
lowStart = datenum('Jul 15, 2009');
lowEnd = lowStart + 340;
lowTicks = datenum(['Aug 1, 2009'; 'Oct 1, 2009'; 'Dec 1, 2009'; ...
    'Feb 1, 2010'; 'Apr 1, 2010'; 'Jun 1, 2010']);
highStart = datenum('Jul 15, 2004');
highEnd = highStart + 340;
highTicks = datenum(['Aug 1, 2004'; 'Oct 1, 2004'; 'Dec 1, 2004'; ...
    'Feb 1, 2005'; 'Apr 1, 2005'; 'Jun 1, 2005']);
testLow_d = decday_d > lowStart & decday_d < lowEnd;
testHigh_d = decday_d > highStart & decday_d < highEnd;
testLow_h = decday_h > lowStart & decday_h < lowEnd;
testHigh_h = decday_h > highStart & decday_h < highEnd;

figure1 = figure(1);
set(figure1, 'Name', ['Site ' num2str(siteID) ' - ' ...
    ' 2 year comparison']);
% Ts and SWE
subplot(2,2,1);
plot(decday_d(testHigh_d), wteq(testHigh_d), 'b', 'LineWidth', 1.5);
xlim([highStart highEnd]); ylim([0 700]);
ylabel('SWE (mm)');
set(gca,'XTick', highTicks, 'XTickLabel', '',...
    'Position',[0.10 0.6 0.39 0.25]);
title('2004-5');
%
subplot(2,2,2);
plot(decday_d(testLow_d), wteq(testLow_d), 'b');
xlim([lowStart lowEnd]); ylim([0 700]);
ylabel('SWE (mm)');
set(gca,'XTick', lowTicks, 'XTickLabel', '', 'YTickLabel', '', ...
    'Position',[0.52 0.6 0.39 0.25]);
title('2009-10');
%
subplot(2,2,3);
plot(decday_d(testHigh_d), airT(testHigh_d), 'Color', [0.7,0.7,0.7],...
    'LineWidth', 1.5);
hold on;
plot(decday_h(testHigh_h), ts(testHigh_h), 'r');
ylim([-25 22]); xlim([highStart highEnd]);
% for some reason the axis changes randomly in the following lines so
% be sure to set the xlimits (above)
zeroline = line(get(gca, 'XLim'), [0, 0]);
set(zeroline, 'Color', 'k', 'LineStyle', ':');
legend('Air Temp', '24hr Mean Ts' );
ylabel('Temp (^oC)');
set(gca,'XTick', highTicks,'Position',[0.10 0.20 0.39 0.4]);
datetick('x','mmm dd', 'keeplimits', 'keepticks');

subplot(2,2,4);
plot(decday_d(testLow_d), airT(testLow_d), 'Color', [0.7,0.7,0.7],...
    'LineWidth', 1.5);
hold on;
plot(decday_h(testLow_h), ts(testLow_h), 'r');
ylim([-25 22]); xlim([lowStart lowEnd]);
zeroline = plot(get(gca, 'xlim'), [0, 0]);
set(zeroline, 'Color', 'k', 'LineStyle', ':');
ylabel('Temp (^oC)');
set(gca,'XTick', lowTicks, 'YTickLabel', '',...
    'Position', [0.52 0.20 0.39 0.4]);
datetick('x', 'mmm dd', 'keeplimits', 'keepticks');

clear all;

% -----------------------------------------------------------------------
% Examine the interannual variability in Ts and vwc at one site
% Load hourly and daily data
siteID = 828; % 828=TrialLake, 972=LouisMeadow, 432=CurrantCreek
              % 330=BeaverDivide, 333=BenLomTrail, 674=OrchardRangeID
              % 654=MudFlatID, 310=BaldyAZ, 720=RockCreek
hourlyData = loadsnotel(siteID, 'hourly');
dailyData = loadsnotel(siteID, 'daily', 'exclude');

% Parse out the date/times
decday_h = datenum(strcat(hourlyData{2}, hourlyData{3}), 'yyyy-mm-ddHH:MM');
wyears_h = unique(hourlyData{10});
startdays = datenum(wyears_h(:)-1, 9, 30);
decday_d = datenum(dailyData{2}, 'yyyy-mm-dd');

% Parse out some variables
wteq = dailyData{4} * 25.4;
%airT = dailyData{9};
ts = filterseries(hourlyData{8}, 'sigma', 25, 3); %{7}=5cm, {8}=20cm, etc
vwc = filterseries(hourlyData{5}, 'sigma', 25, 3);

% Overlap all years of Tair and Tsoil data (AGU 2011 poster)
h = figure(2);
set(h, 'Name', ['Site ' num2str(siteID) ' - ' ...
    ' 2 year comparison']);

% X Tick locations and labels
ticklocs = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365];
tickmonths = ['Oct'; 'Nov'; 'Dec'; 'Jan'; 'Feb'; 'Mar'; 'Apr'; 'May';...
    'Jun';'Jul'; 'Aug'; 'Sep'; 'Oct'];

subplot(3,1,1);
hold on;
% Initialize variables for calculating a mean timeseries
wteqConcat = []; doyConcat = [];
for i = 1:length(wyears_h)
    % Slice out each wateryear worth of data/datenums
    dailytest = dailyData{21}==wyears_h(i);
    wyWteq = wteq(dailytest);
    % Subtract of the initial datenum to get doy
    doys = decday_d(dailytest) - startdays(i);
    plot(doys, wyWteq, 'Color', [0.7,0.7,0.7], ...
        'LineWidth', 1.5);
    wteqConcat = [wteqConcat; wyWteq]; % Concatenate yearly wteq
    doyConcat = [doyConcat; doys]; % And each years doy values
end;
doyConcat = [doyConcat ones(size(doyConcat))]; % Create accumarray index
% Get a mean timeseries with accumarray and plot it
wteqMean = accumarray(doyConcat, wteqConcat, ...
    [numel(unique(doyConcat)) 1], @nanmean)
plot(1:366, wteqMean, '-b', 'LineWidth', 2);
% Set axes limits, tick locations, labels, position, etc
xlim([0 367]); ylim([-5 1500]);
ylabel('mm');
set(gca,'XTick',ticklocs, 'XTickLabel', '',...
    'Position', [0.13, 0.678, 0.775, 0.25]);
title('Trial Lake');
%
subplot(3,1,2);
hold on;
tsConcat = []; doyConcat = [];
for i = 1:length(wyears_h)
    hourlytest = hourlyData{10}==wyears_h(i);
    wyTs = ts(hourlytest);
    doys = decday_h(hourlytest) - startdays(i);
    plot(doys, wyTs, 'Color', [0.7,0.7,0.7],...
        'LineWidth', 1.5);
    tsConcat = [tsConcat; wyTs];
    doyConcat = [doyConcat; doys];
end;
[doyvals, ~, doyindex] = unique(doyConcat);
doyConcat = [doyindex ones(size(doyindex))];
tsMean = accumarray(doyConcat, tsConcat,...
    [numel(unique(doyConcat)) 1], @nanmean);
plot(doyvals, tsMean, '-r');
% Set axes limits, tick locations, labels, position, etc
zeroline = line(get(gca, 'XLim'), [0, 0]);
set(zeroline, 'Color', 'k', 'LineStyle', ':');
xlim([0 367]); ylim([-5 25]);
ylabel('^oC');
set(gca,'XTick',ticklocs, 'XTickLabel', '',...
    'Position', [0.13, 0.424, 0.775, 0.25]);

subplot(3,1,3);
hold on;
vwcConcat = []; doyConcat = [];
for i = 1:length(wyears_h)
    hourlytest = hourlyData{10}==wyears_h(i);
    wyVwc = vwc(hourlytest);
    doys = decday_h(hourlytest) - startdays(i);
    plot(doys, vwc(hourlytest), 'Color', [0.7,0.7,0.7],...
        'LineWidth', 1.5);
    vwcConcat = [vwcConcat; wyVwc];
    doyConcat = [doyConcat; doys];
end;
[doyvals, ~, doyindex] = unique(doyConcat);
doyConcat = [doyindex ones(size(doyindex))];
vwcMean = accumarray(doyConcat, vwcConcat,...
    [numel(unique(doyConcat)) 1], @nanmean);
plot(doyvals, vwcMean, '-k', 'LineWidth', 2);
% Set axes limits, tick locations, labels, position, etc
xlim([0 367]); ylim([-5 45]);
ylabel('VWC (%)');
set(gca,'XTick', ticklocs, 'XTickLabel', tickmonths,...
    'Position', [0.13, 0.172, 0.775, 0.25]);
    
