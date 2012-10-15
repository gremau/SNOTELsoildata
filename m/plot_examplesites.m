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


% -----------------------------------------------------------------------
% Examine the seasonal variability in temp or vwc at a set of 4 sites.
% These sites have been chosen to represent elevation/temp and SWE
% gradients
%
% 2 plots = Sensor timeseries for each site and then seasonal histograms

% Set list of sites, sensor output(vwc or temp), and sensor depth
siteIDs = [828, 333, 452, 336]; % elev/swe: hi/hi, low/hi, hi/low, low/low
% 828 = TrialLk, 333 = BenLomTrail, 452=DonkeyRes, 573=Big BendNV
sensoroutput = 'vwc';
sensordepth = 2; %(1=5cm, 2=20cm, 3=50cm);
startwy = 2006;

% Select TEMP or VWC data and set distribution bins and plot axes
if strcmpi(sensoroutput, 'vwc');
    sensorcolumn = sensordepth + 3; % get proper column using sensordepth
    xmin = 0;
    % If running RAW SENSOR DATA (no normalization)
    % xedges = 0:1:100; % raw sm data bins (0-100)
    % xmax = 75
    % If running NORMALIZED data with smnormalize
    xedges = 0:0.01:1; % normalized vwc bins (0-1)
    xmax = 1; % these axes are good for normalized data
    ymax = 0.125;
    disp('*** Running in normalized soil moisture data mode ***');
elseif strcmpi(sensoroutput, 'temp');
    sensorcolumn = sensordepth + 6; % get proper column using sensordepth
    xedges = -15:0.5:35; % soil temp bins -15-35 deg C
    xmax = 35; % corresponding x axis limits
    xmin = -15;
    ymax = 0.2;
end

% Set up PLOT 1 - add each site's timeseries on iteration through following
% loop
h = figure(3);

% Allocate for histogram and mean matrices - fill in following loop
histograms = zeros(length(xedges), 16);
means = zeros(16, 1);

for i = 1:length(siteIDs);
    % Load hourly data from site  w/ loadsnotel:
    siteHourly = loadsnotel(siteIDs(i), 'hourly', 'exclude');
    % Get rid of wateryears prior to startwy
    wyexclude = siteHourly{10}>startwy-1;
    for j = 1:10
        siteHourly{j} = siteHourly{j}(wyexclude);
    end
    % Parse out the desired sensor depth, normalize if plotting vwc
    if strcmpi(sensoroutput, 'vwc')
        sensordata = filterseries(siteHourly{sensorcolumn}, 'sigma', 25, 3);
        % SPECIAL Case for Taylor Cyn - There is some bad data that makes it
        % past filter and messes up normalization - remove it
        if siteIDs(i) == 336 %811n for TaylorCyn, 336 for BigBend
            test = sensordata<10; %18 for TaylorCyn, 10 for BigBend
            sensordata(test)=nan;
        end
        sensordata = smnormalize(sensordata, 1);
    else
        sensordata = siteHourly{sensorcolumn};
    end
    
    % Create date arrays
    datevec_h = datevec(strcat(siteHourly{2}, siteHourly{3}), 'yyyy-mm-ddHH:MM');
    datenum_h = datenum(datevec_h);
    
    % PLOT 1 - add the entire timeseries for site i
    ticklocations = linspace(min(datenum_h), max(datenum_h), 20);
    set(h, 'Name', ['Site ' num2str(siteIDs(i)) ' - Full sensor timeseries']);
    subplot (4, 1, i)
    plot(datenum_h, sensordata, 'k');
    set(gca, 'XTick', ticklocations);
    set(gca, 'XTickLabel', ticklocations);
    datetick('x', 12, 'keepticks');
    ylabel(sensoroutput);
    title(['Site ' num2str(siteIDs(i))]);
    
    % Create logical tests and pull desired quarters (3 months intervals)
    testOND = (datevec_h(:,2)==10 | datevec_h(:,2)==11 | datevec_h(:,2)==12);
    testJFM = (datevec_h(:,2)==1 | datevec_h(:,2)==2 | datevec_h(:,2)==3);
    testAMJ = (datevec_h(:,2)==4 | datevec_h(:,2)==5 | datevec_h(:,2)==6);
    testJAS = (datevec_h(:,2)==7 | datevec_h(:,2)==8 | datevec_h(:,2)==9);
    sensordata_OND = sensordata(testOND,:);
    sensordata_JFM = sensordata(testJFM,:);
    sensordata_AMJ = sensordata(testAMJ,:);
    sensordata_JAS = sensordata(testJAS,:);
    
    % Generate histograms for each quarter's sensor data
    histOND = histc(sensordata_OND, xedges);
    histJFM = histc(sensordata_JFM, xedges);
    histAMJ = histc(sensordata_AMJ, xedges);
    histJAS = histc(sensordata_JAS, xedges);
    
    % Normalize and put histograms in the histograms matrix
    % (in plotting order)
    histograms(:, i) = histOND./sum(histOND);
    histograms(:, i+4) = histJFM./sum(histJFM);
    histograms(:, i+8) = histAMJ./sum(histAMJ);
    histograms(:, i+12) = histJAS./sum(histJAS);
    
    % Calculate mean and standard deviations of data from each quarter and
    % place in appropriate vector (in plotting order)
    means(i) = mean(sensordata_OND(~isnan(sensordata_OND)));
    means(i+4) = mean(sensordata_JFM(~isnan(sensordata_JFM)));
    means(i+8) = mean(sensordata_AMJ(~isnan(sensordata_AMJ)));
    means(i+12) = mean(sensordata_JAS(~isnan(sensordata_JAS)));
    % stddevOND = std(sensordata_OND(~isnan(sensordata_OND)));
    % stddevJFM = std(sensordata_JFM(~isnan(sensordata_JFM)));
    % stddevAMJ = std(sensordata_AMJ(~isnan(sensordata_AMJ)));
    % stddevJAS = std(sensordata_JAS(~isnan(sensordata_JAS)));

end
    clear testOND testJFM testAMJ testJAS;

% PLOT 2. Plot quarterly distributions for all sites
titles = {['Site ' num2str(siteIDs(1)) ' Oct-Dec'] ...
    ['Site ' num2str(siteIDs(2)) ' Oct-Dec']...
    ['Site ' num2str(siteIDs(3)) ' Oct-Dec']...
    ['Site ' num2str(siteIDs(4)) ' Oct-Dec'] ...
    'Jan-Mar' 'Jan-Mar' 'Jan-Mar' 'Jan-Mar' ...
    'Apr-Jun' 'Apr-Jun' 'Apr-Jun' 'Apr-Jun' ...
    'Jul-Sep' 'Jul-Sep' 'Jul-Sep' 'Jul-Sep'};

h = figure(4);
set(h, 'Name', ['4 SNOTEL Sites - ' sensoroutput ...
    ' quarterly histograms - all years combined']);
% Loop through 16 subplots and plot histograms and means
for i = 1:16;
    subplot (4, 4, i)
    bar (xedges, histograms(:, i), 'g');
    hold on
    plot([means(i) means(i)], [0 1], ':k');
    axis([xmin xmax 0 ymax]);
    ylabel('Frequency');
    title(titles{i});
end
    
