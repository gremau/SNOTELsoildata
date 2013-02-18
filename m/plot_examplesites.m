% plot_examplesites.m
%
% This makes a number of figures that highlight snow-soil interactions
% at particular sites.
% **File input:** The first 4 figures directly load hourly or daily data. 
%                 The last 2 use the summary files
% **User input:** User can select different sites to plot in figs 2, 4, 
%                 5, 6 in the code
% 
% Fig 1 = 2 year SWE/temp comparison at Mosby Mtn.
% Fig 2 = Multi-year SWE/Tsoil/VWC time series at a site (changeable)
% Fig 3 = Full VWC timeseries of 4 sites - the data for Fig 4
% Fig 4 = VWC Frequency histograms for 4 contrasting sites
% Fig 5 = Scatterplot/Regression of snowcov. Tsoil vs pre-onset T at 1 site
% Fig 6 = Scatterplot/Regression of summer (JAS) VWC vs meltday at 1 site

close all;
clear all;
addpath('/home/greg/data/code_resources/m_common/linreg/'); 
addpath('/home/greg/data/code_resources/m_common/');


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
sweyticks = [200; 400; 600];
testLow_d = decday_d > lowStart & decday_d < lowEnd;
testHigh_d = decday_d > highStart & decday_d < highEnd;
testLow_h = decday_h > lowStart & decday_h < lowEnd;
testHigh_h = decday_h > highStart & decday_h < highEnd;

figure1 = figure('position',[100 0 1100 800],'paperpositionmode',...
    'auto', 'color','none','InvertHardcopy','off');

set(figure1, 'Name', ['Site ' num2str(siteID) ' - ' ...
    ' 2 year comparison']);
set(figure1, 'DefaultAxesFontSize',16, 'DefaultTextFontSize', 18);

% Ts and SWE
subplot(2,2,1);
plot(decday_d(testHigh_d), wteq(testHigh_d), '-k', 'LineWidth', 2);
%plot(decday_d(testHigh_d), wteq(testHigh_d), 'Color', rgb('Gold'),...
 %   'LineWidth', 2);
xlim([highStart highEnd]); ylim([0 700]);
ylabel('SWE (mm)');
set(gca,'XTick', highTicks, 'XTickLabel', '','YTick', sweyticks,...
    'Position',[0.10 0.6 0.39 0.25]);
title('2004-5', 'Fontsize', 20, 'Fontangle', 'italic');
%
subplot(2,2,2);
plot(decday_d(testLow_d), wteq(testLow_d), '-k', 'LineWidth', 2);
%plot(decday_d(testLow_d), wteq(testLow_d), 'Color', rgb('Gold'),...
 %   'LineWidth', 2);
xlim([lowStart lowEnd]); ylim([0 700]);
set(gca,'XTick', lowTicks, 'XTickLabel', '', 'YTickLabel', '', ...
    'Position',[0.52 0.6 0.39 0.25]);
title('2009-10', 'Fontsize', 20, 'Fontangle', 'italic');
%
subplot(2,2,3);
plot(decday_d(testHigh_d), airT(testHigh_d), 'Color', [0.5,0.5,0.5],...
    'LineWidth', 1.5);
hold on;
plot(decday_h(testHigh_h), ts(testHigh_h), '-k', 'LineWidth', 1.5)
%plot(decday_h(testHigh_h), ts(testHigh_h), 'Color', rgb('DarkRed'),...
 %'LineWidth', 1.5);
ylim([-22 22]); xlim([highStart highEnd]);
% for some reason the axis changes randomly in the following lines so
% be sure to set the xlimits (above)
zeroline = line(get(gca, 'XLim'), [0, 0]);
set(zeroline, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1.5);
legend('T_{air}', 'T_{soil}', 'location', 'southwest' );
ylabel('T (^oC)');
set(gca,'XTick', highTicks,'Position',[0.10 0.20 0.39 0.4]);
datetick('x','mmm', 'keeplimits', 'keepticks');

subplot(2,2,4);
plot(decday_d(testLow_d), airT(testLow_d), 'Color', [0.5,0.5,0.5],...
    'LineWidth', 1.5);
hold on;
plot(decday_h(testLow_h), ts(testLow_h), '-k',...
    'LineWidth', 1.5);
%plot(decday_h(testLow_h), ts(testLow_h), 'Color', rgb('DarkRed'),...
 %   'LineWidth', 1.5);
ylim([-22 22]); xlim([lowStart lowEnd]);
zeroline = plot(get(gca, 'xlim'), [0, 0]);
set(zeroline, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1.5);
set(gca,'XTick', lowTicks, 'YTickLabel', '',...
    'Position', [0.52 0.20 0.39 0.4]);
datetick('x', 'mmm', 'keeplimits', 'keepticks');

figpath = '../figures/';
print(figure1,'-depsc2','-painters',[figpath 'figF.eps']) 
clear all;

% -----------------------------------------------------------------------
% Examine the interannual variability in Ts and vwc at one site
% Load hourly and daily data
siteID = 432; % 828=TrialLake, 972=LouisMeadow, 432=CurrantCreek
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
Tair = dailyData{9};
ts = filterseries(hourlyData{8}, 'sigma', 25, 3); %{7}=5cm, {8}=20cm, etc
vwc = filterseries(hourlyData{5}, 'sigma', 25, 3);

% Overlap all years of Tair and Tsoil data (AGU 2011 poster)
figure2 = figure('position',[100 0 1100 800],'paperpositionmode',...
    'auto', 'color','none','InvertHardcopy','off');
set(figure2, 'Name', ['Site ' num2str(siteID) ' - ' ...
    ' 2 year comparison']);
set(figure2, 'DefaultAxesFontSize',16, 'DefaultTextFontSize', 18);

% X Tick locations and labels
ticklocs = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365];
tickmonths = ['Oct'; 'Nov'; 'Dec'; 'Jan'; 'Feb'; 'Mar'; 'Apr'; 'May';...
    'Jun';'Jul'; 'Aug'; 'Sep'; 'Oct'];

subplot(3,1,1);
% Initialize variables for calculating a mean timeseries
wteqConcat = []; doyConcat = [];
for i = 1:length(wyears_h)
    % Slice out each wateryear worth of data/datenums
    dailytest = dailyData{21}==wyears_h(i);
    wyWteq = wteq(dailytest);
    % Subtract of the initial datenum to get doy
    doys = decday_d(dailytest) - startdays(i);
    plot(doys, wyWteq, 'Color', [0.5,0.5,0.5], 'Linewidth', 1.5);
    hold on;
    wteqConcat = [wteqConcat; wyWteq]; % Concatenate yearly wteq
    doyConcat = [doyConcat; doys]; % And each years doy values
end;
doyConcat = [doyConcat ones(size(doyConcat))]; % Create accumarray index
% Get a mean timeseries with accumarray and plot it
wteqMean = accumarray(doyConcat, wteqConcat, ...
    [numel(unique(doyConcat)) 1], @nanmean);
plot(1:366, wteqMean, '-k', 'LineWidth', 2);
%plot(1:366, wteqMean, '-', 'Color', rgb('Gold'),'LineWidth', 2);
% Set axes limits, tick locations, labels, position, etc
xlim([0 367]); ylim([-5 450]);
ylabel('SWE (mm)');
set(gca,'XTick',ticklocs, 'XTickLabel', '', 'Ytick', [100;200;300;400],...
    'Position', get(gca, 'position') .* [1 .9 1 1.23]);
text(0.95, 0.85, 'A', 'Units', 'normalized');
%title('Currant Creek, UT', 'Fontsize', 20, 'Fontangle', 'italic');
%
subplot(3,1,2);
tsConcat = []; doyConcat = [];taConcat = []; doy_dConcat = [];
for i = 1:length(wyears_h)
    hourlytest = hourlyData{10}==wyears_h(i);
    wyTs = ts(hourlytest);
    doys = decday_h(hourlytest) - startdays(i);
    plot(doys, wyTs, 'Color', [0.5,0.5,0.5], 'Linewidth', 1.5);
    hold on;
    tsConcat = [tsConcat; wyTs];
    doyConcat = [doyConcat; doys];
    %Do this for Tair also
    dailytest = dailyData{21}==wyears_h(i);
    wyTair = Tair(dailytest);
    doys_d = decday_d(dailytest) - startdays(i);
    taConcat = [taConcat; wyTair];
    doy_dConcat = [doy_dConcat; doys_d];
end;
% Plot Tsoil mean
[doyvals, ~, doyindex] = unique(doyConcat);
doyConcat = [doyindex ones(size(doyindex))];
tsMean = accumarray(doyConcat, tsConcat,...
    [numel(unique(doyConcat)) 1], @nanmean);
h1 = plot(doyvals, tsMean, '-k', 'Linewidth', 1.5);
%h1 = plot(doyvals, tsMean, '-', 'Color', rgb('DarkRed'), 'Linewidth', 2);
% Plot Tair mean
[doyvals, ~, doyindex] = unique(doy_dConcat);
doy_dConcat = [doyindex ones(size(doyindex))];
taMean = accumarray(doy_dConcat, taConcat,...
    [numel(unique(doy_dConcat)) 1], @nanmean);
%h2 = plot(doyvals, taMean, '--k', 'Linewidth', 2);
% Set axes limits, tick locations, labels, position, etc
zeroline = line(get(gca, 'XLim'), [0, 0]);
set(zeroline, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1.5);
xlim([0 367]); ylim([-5 20]);
ylabel('T_{soil} (^oC)');
%legend([h1 h2], 'T_{soil}', 'T_{air}', 'location', 'southeast');
set(gca,'XTick',ticklocs, 'XTickLabel', '',...
    'Position', get(gca, 'position') .* [1 .9 1 1.23]);
text(0.95, 0.85, 'B', 'Units', 'normalized');

subplot(3,1,3);
vwcConcat = []; doyConcat = [];
for i = 1:length(wyears_h)
    hourlytest = hourlyData{10}==wyears_h(i);
    wyVwc = vwc(hourlytest);
    doys = decday_h(hourlytest) - startdays(i);
    plot(doys, vwc(hourlytest), 'Color', [0.5,0.5,0.5], 'Linewidth', 1.5);
    hold on;
    vwcConcat = [vwcConcat; wyVwc];
    doyConcat = [doyConcat; doys];
end;
[doyvals, ~, doyindex] = unique(doyConcat);
doyConcat = [doyindex ones(size(doyindex))];
vwcMean = accumarray(doyConcat, vwcConcat,...
    [numel(unique(doyConcat)) 1], @nanmean);
%plot(doyvals, vwcMean, '-', 'Color', rgb('Navy'),'LineWidth', 2);
plot(doyvals, vwcMean, '-k','LineWidth', 1.5);
zeroline = line(get(gca, 'XLim'), [0, 0]);
set(zeroline, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1.5);
% Set axes limits, tick locations, labels, position, etc
xlim([0 367]); ylim([-2 35]);
ylabel('VWC (%)');
set(gca,'XTick', ticklocs, 'XTickLabel', tickmonths,...
    'Position', get(gca, 'position') .* [1 .9 1 1.23]);
text(0.95, 0.85, 'C', 'Units', 'normalized');

figpath = '../figures/';
print(figure2,'-depsc2','-painters',[figpath 'figB.eps'])

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
    xedges = 0:0.02:1; % normalized vwc bins (0-1)
    xmax = 1; % these axes are good for normalized data
    ymax = 0.25;
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
figure3 = figure();

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
    datevec_h = datevec(strcat(siteHourly{2}, siteHourly{3}),...
        'yyyy-mm-ddHH:MM');
    datenum_h = datenum(datevec_h);
    
    % PLOT 1 - add the entire timeseries for site i
    ticklocations = linspace(min(datenum_h), max(datenum_h), 20);
    set(figure3, 'Name', ['Site ' num2str(siteIDs(i))...
        ' - Full sensor timeseries']);
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
titlelabels = {'Hi SWE/Hi Elev' 'Hi SWE/Low Elev' 'Low SWE/Hi Elev'...
    'Low SWE/Low Elev'};
monthlabels = {'Oct-Dec' '' '' '' 'Jan-Mar' '' '' '' 'Apr-Jun' '' '' ''...
    'Jul-Sep'};

figure4 = figure('position',[100 0 1000 800],'paperpositionmode',...
    'auto', 'color','none','InvertHardcopy','off');
set(figure4, 'Name', ['4 SNOTEL Sites - ' sensoroutput ...
    ' quarterly histograms - all years combined']);
set(figure4, 'DefaultAxesFontSize',16, 'DefaultTextFontSize', 16);

% Loop through 16 subplots and plot histograms and means
for i = 1:16;
    subplot (4, 4, i)
    bar (xedges, histograms(:, i), 'Facecolor', [0.7 0.7 0.7]);
    hold on
    plot([means(i) means(i)], [0 1], '--k', 'Linewidth', 1.5);
    axis([xmin xmax 0 ymax]);
    set(gca, 'position', [0.925 0.925 1.15 1.19] .* get(gca, 'position'),...
        'Xtick', [0.5]);
    if i==1 || i==5 || i==9 || i==13;
        %ylabel('Frequency');
        text(0.1, 0.8, monthlabels(i), 'Units', 'normalized');
        set(gca, 'XTick', [0;0.5])
        if i==5;
            text(-0.35, -1,'Normalized frequency of ocurrence',...
                'Units', 'normalized', 'Rotation', 90);
        end
    elseif i==16;
        set(gca, 'Xtick', [0.5;1]);
        set(gca, 'YtickLabel', '');
    else
        set(gca, 'YtickLabel', '');
    end
    if i < 5
        %title({['Site ' num2str(siteIDs(i))]; titlelabels{i}},...
        %    'Fontsize', 18, 'Fontangle', 'italic');
        title(titlelabels{i}, 'Fontsize', 18, 'Fontangle', 'italic');
        set(gca, 'XtickLabel', '');
    elseif i < 13
        set(gca, 'XtickLabel', '');
    elseif i 
    end
end

figpath = '../figures/';
print(figure4,'-depsc2','-painters',[figpath 'figJ.eps'])

% -----------------------------------------------------------------------
% Look at linear regressions of Mean annual soil temperature and mean
% summer soil moisture versus snowpack at 2 sites. These are examples to
% put regression results in context.
%
% 2 plots = One site regression for MAST and one for JAS VWC

% Add any needed tools
addpath('/home/greg/data/code_resources/m_common/'); 
addpath('/home/greg/data/code_resources/m_common/nanstuff/');
addpath('/home/greg/data/code_resources/m_common/linear/'); 
addpath('/home/greg/data/code_resources/m_common/hline_vline/'); 

% Use this site as an example:
% mastSite = 828;
% mastSitelabel = 'Trial Lake';
mastSite = 393;
mastSitelabel = 'Chalk Creek';
% mastSite = 332;
% mastSitelabel = 'Ben Lomond Peak';
% vwcSite = 333;
% vwcSitelabel = 'Ben Lomond Trail';
% vwcSites = 390;
% vwcSitelabel = 'Castle Valley';
% vwcSite = 766;
% vwcSitelabel = 'Snowbird';
vwcSite = 839;
vwcSitelabel = 'Upper Rio Grande';

% Set processed data path
processeddatapath = '../processed_data/';

% LOAD the data (can switch between daily/hourly data here)
climData = csvread([processeddatapath 'wyear_climatesummary.txt'],1,0);
tsData = csvread([processeddatapath 'wyear_soiltempsummary_hourly.txt'],1,0);
% tsData = csvread([processeddatapath 'wyear_soiltempsummary_daily.txt'],1,0);
vwcDataN = csvread([processeddatapath ...
    'wyear_soilwatersummary_hourly_smnorm.txt'],1,0);

% Get a subset of climData that corresponds with available soildata
[matchsoil, idx] = ismember(climData(:, 1:2), tsData(:, 1:2), 'rows');
soilClim = climData(matchsoil, :);
% matchsoil2 = ismember(tsData(:, 1:2), soilClim(:, 1:2), 'rows');

% Climate data
site_cl = soilClim(:, 1);
year_cl = soilClim(:, 2);
maxswe = soilClim(:, 3)*25.4;
maxsweday = soilClim(:, 4);
maxdepth = soilClim(:, 5);
onsetdoy = soilClim(:, 6);
meltdoy = soilClim(:, 7);
snowduration = soilClim(:, 8);
totaldaysSC = soilClim(:, 9); % Total days, may vary from duration above
maxcontinSC = soilClim(:, 10);% Length of longest continuos snowpack
numcontinSC = soilClim(:, 11);% # of continuous snowcovered periods
accumprecip = soilClim(:, 12)*25.4;
JASprecip = soilClim(:, 13)*25.4;
maat = soilClim(:, 80);
maat_sd = soilClim(:, 81);

preonsetTair = soilClim(:, 82);
preonsetTairSd = soilClim(:, 83);
premeltTair = soilClim(:, 84);
premeltTairSd = soilClim(:, 85);
postmeltTair = soilClim(:, 86);
postmeltTairSd = soilClim(:, 87);
elev = soilClim(:, 88);
lat = soilClim(:, 89);
lon = soilClim(:, 90);
ltMeanSWE = soilClim(:, 91);
ltMeanPrecip = soilClim(:, 92);


% Seasonal/yearly soil temp means
site_ts = tsData(:, 1);
mast5cm = tsData(:, 99);
sdast5cm = tsData(:, 100);
mast20cm = tsData(:, 101);
sdast20cm = tsData(:, 102);
mast50cm = tsData(:, 103);
sdast50cm = tsData(:, 104);
% Snowcovered soil temp means
snowcovTs5mean = tsData(:, 105);
snowcovTs5sd = tsData(:, 106);
snowcovTs20mean = tsData(:, 107);
snowcovTs20sd = tsData(:, 108);
snowcovTs50mean = tsData(:, 109);
snowcovTs50sd = tsData(:, 110);

% Seasonal soil moisture
amjVWC5mean = vwcDataN(:, 87);
amjVWC5sd = vwcDataN(:, 88);
amjVWC20mean = vwcDataN(:, 89);
amjVWC20sd = vwcDataN(:, 90);
amjVWC50mean = vwcDataN(:, 91);
amjVWC50sd = vwcDataN(:, 92);
jasVWC5mean = vwcDataN(:, 93);
jasVWC5sd = vwcDataN(:, 94);
jasVWC20mean = vwcDataN(:, 95);
jasVWC20sd = vwcDataN(:, 96);
jasVWC50mean = vwcDataN(:, 97);
jasVWC50sd = vwcDataN(:, 98);

preonsetVWC5 = vwcDataN(:, 99);
preonsetVWC5sd = vwcDataN(:, 100);
preonsetVWC20 = vwcDataN(:, 101);
preonsetVWC20sd = vwcDataN(:, 102);
preonsetVWC50 = vwcDataN(:, 103);
preonsetVWC50sd = vwcDataN(:, 104);

sites = unique(site_cl);

% The variables to use  
diff5cm = mast5cm-maat; diff20cm = mast20cm-maat; diff50cm = mast50cm-maat;
inter = preonsetTair .* onsetdoy;
Yvars = {'snowcovTs20mean' 'jasVWC20mean'};
Xvars = {'inter' 'meltdoy'};


% FIG 1 - Plot the results for the MAST site
figure5 = figure('position',[100 0 600 400],'paperpositionmode',...
    'auto', 'color','none','InvertHardcopy','off');
set(figure5, 'DefaultAxesFontSize',16, 'DefaultTextFontSize', 16);
set(figure5, 'Name', ['Regression: Y = ' Yvars{1} ', X = ' Xvars{1} ...
    ', Site = ' num2str(mastSite)]);
test = site_cl==mastSite;
% Be sure that the variables here are used for all sites above
eval(['ysite = ' Yvars{1} '(test);']);
eval(['xsite = ' Xvars{1} '(test);']);

plot(xsite, ysite, 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'Black');
hold on;
xrange = xlim(gca);
[coeffs, rsq, xfit, yfit] = fitline(xsite, ysite, 1, xrange);
[b,bint,resid,rint,stats] = regress2(ysite, [xsite ones(size(xsite))]);
plot(xfit, yfit,'--k', 'LineWidth', 1.5);
text(0.6, 0.8,['r^2 = ' num2str(rsq, 2) ', p = ' num2str(stats(3), 2)],...
    'Units', 'Normalized'); % r^2 & p
xlabel('Number of snow-covered days'); ylabel('Mean annual T_{soil} (^oC)');
title(mastSitelabel);

figpath = '../figures/';
print(figure5,'-depsc2','-painters',[figpath 'figE.eps'])
    
% FIG 2 - Plot the results for the VWC site
figure6 = figure('position',[100 0 600 400],'paperpositionmode',...
    'auto', 'color','none','InvertHardcopy','off');
set(figure6, 'DefaultAxesFontSize',16, 'DefaultTextFontSize', 16);
set(figure6, 'Name', ['Regression: Y = ' Yvars{2} ', X = ' Xvars{2} ...
    ', Site = ' num2str(vwcSite)]);
test = site_cl==vwcSite;
% Be sure that the variables here are used for all sites above
eval(['ysite = ' Yvars{2} '(test);']);
eval(['xsite = ' Xvars{2} '(test);']);

plot(xsite, ysite, 'ok', 'MarkerSize', 10, 'MarkerFaceColor', 'Black');
hold on;
xrange = xlim(gca);
[coeffs, rsq, xfit, yfit] = fitline(xsite, ysite, 1, xrange);
[b,bint,resid,rint,stats] = regress2(ysite, [xsite ones(size(xsite))]);
plot(xfit, yfit,'--k', 'LineWidth', 1.5);
text(0.6, 0.8,['r^2 = ' num2str(rsq, 2) ', p = ' num2str(stats(3), 2)],...
    'Units', 'Normalized'); % r^2 & p
xlabel('Day of Snowmelt'); ylabel('Mean summer VWC (20cm)');
title(vwcSitelabel);

figpath = '../figures/';
print(figure6,'-depsc2','-painters',[figpath 'figH.eps'])
