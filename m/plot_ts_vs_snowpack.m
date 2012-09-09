% plot_ts_vs_snowpack.m
%
% Plots soil and air temp from recent daily sensor data vs SWE, snow 
% depth for SNOTEL sites (and subsets of SNOTELS)
%
% Uses daily sensor data in the rawdata/SNOTEL data/ directory
%
% Version 1: 111116
%

clear;          % clear memory
close all;      % clear any figures
fignum=0;       % used to increment figure number for plots
%addpath('../m/');
addpath('~/data/code_resources/m_common/stat_tbx/');

% Ask user for month number
monthsel = str2double(input('Which month (1-12)?: ', 's'));

% Set data path and file name, read in file
rawdatapath = '../rawdata/';
processeddatapath = '../processed_data/';

% Load list of sites with data in the daily data directory
dailyDataSites = sortrows(csvread([rawdatapath 'allsensors_daily/sitelist.txt']));

% Import list of wasatch + uinta sites
formatstr = '%s%f%s%s';
fid = fopen([processeddatapath 'SNOTELrangelist.csv']);
wasatchUintaCell = textscan(fid, formatstr,'Headerlines', 1, 'Delimiter', ',');
fclose(fid);

% Creat list of wasatch and uinta sites
wasatchTest = strcmpi(wasatchUintaCell{4}, 'WASATCH');
uintaTest =  strcmpi(wasatchUintaCell{4},'UINTA');
wasatch = wasatchUintaCell{2}(wasatchTest);
uintas = wasatchUintaCell{2}(uintaTest);
clear test;

sites = unique(dailyDataSites(:, 1));
monthLabels = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sept' 'Oct'...
    'Nov' 'Dec'};
monthlabel = monthLabels{monthsel};
monthMeans = [];
% Load data and parse out month data
for i = 1:length(sites);
    m = loadsnotel('daily', sites(i));
    %Create datevector for datafile
    siteDateVec = datevec(m{2}, 'yyyy-mm-dd');
    % Columns are site, year, st-5, st-20, st-60, sndepth, swe, airT
    siteData = [double(m{1}) siteDateVec(:,1) m{14} m{15} m{16} (m{10}*25.4) (m{4}*25.4) m{9}];
    % Get monthly data
    monthTest = siteDateVec(:,2)==monthsel;
    monthData = siteData(monthTest, :);
    monthYears = unique(monthData(:,2));
    % Reduce to yearly averages
    for j = 1:length(monthYears)
        yearTest = monthData(:, 2) == monthYears(j);
        monthMeans = [monthMeans; nanmean(monthData(yearTest, :))];
    end
end


% PLOTS
%----------------------------------------------------------------------
% FIG 1 - Month soil temps vs snowpack
% Note that each datapoint is one year of temp and snowpack data at one
% site during the month of interest
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name', ['Mean ' monthlabel ' Ts vs snowpack at 2 depths']);

% Mean month soil temp by SWE - 5cm
subplot 221;
plot(monthMeans(:,7), monthMeans(:,3), 'ok');
xlabel('Mean SWE');
ylabel('Mean 5cm soil temp (Celsius)');
title([monthlabel ' 5cm soil temp vs ' monthlabel ' SWE']);

% Mean month soil temp by snow depth - 5cm
subplot 222;
plot(monthMeans(:,6), monthMeans(:,3), 'ok');
xlabel('Mean snow depth');
%ylabel('Mean soil temp (Celsius)');
title(['vs ' monthlabel ' snow depth']);

% Mean month soil temp by SWE - 20cm
subplot 223;
plot(monthMeans(:,7), monthMeans(:,4), 'ok');
xlabel('Mean SWE');
ylabel('Mean 20cm soil temp (Celsius)');
title([monthlabel ' 20cm soil temp vs ' monthlabel ' SWE']);

% Mean month soil temp by snow depth - 20cm
subplot 224;
plot(monthMeans(:,6), monthMeans(:,4), 'ok');
xlabel('Mean snow depth');
%ylabel('Mean soil temp (Celsius)');
title(['vs ' monthlabel ' snow depth']);

%--------------------------------------------------------
% FIG 2 - Same as above, but tweaked for AGU 2011 poster
fignum = fignum+1;    
h = figure(fignum);

subplot 121;
plot(monthMeans(:,7), monthMeans(:,3), 'ob');
xlabel('Mean SWE (mm)');
ylabel('Soil T (^oC)');
title('December soil temperatures');
legend('5cm one-month mean');

% Mean month soil temp by SWE - 20cm
subplot 122;
plot(monthMeans(:,7), monthMeans(:,4), 'ok');
xlabel('Mean SWE (mm)');
ylabel('Soil T (^oC)');
legend('20cm one-month mean');
%title([monthlabel ' 20cm temp vs ' monthlabel ' SWE']);

%------------------------------------------------------
% FIG 3 - Same as above, but plot for wasatch and uinta sites 
uintaTest = ismember(monthMeans(:, 1), uintas);
wasatchTest = ismember(monthMeans(:, 1), wasatch);

fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name', ['Mean ' monthlabel ' Ts vs snowpack - Wasatch & Uintas']);
% Mean month soil temp by SWE - 5cm
subplot 221;
plot(monthMeans(wasatchTest,7), monthMeans(wasatchTest,3), 'om');
hold on;
plot(monthMeans(uintaTest,7), monthMeans(uintaTest,3), 'ob');
bothTest = wasatchTest | uintaTest;
x = monthMeans(bothTest, 7);
y = monthMeans(bothTest, 3);
yNoNan = ~isnan(y);
coefficients = polyfit(x(yNoNan), y(yNoNan), 2);
pnFit = polyval(coefficients, (0:1:700));
plot((0:1:700), pnFit,':k');
xlabel('Mean SWE');
ylabel('Mean 5cm soil temp (Celsius)');
legend('Wasatch mtns', 'Uinta mtns');
title([monthlabel ' 5cm soil temp vs ' monthlabel ' SWE']);

% Mean month soil temp by snow depth - 5cm
subplot 222;
plot(monthMeans(wasatchTest,6), monthMeans(wasatchTest,3), 'om');
hold on
plot(monthMeans(uintaTest,6), monthMeans(uintaTest,3), 'ob');
xlabel('Mean snow depth');
%ylabel('Mean soil temp (Celsius)');
title(['vs ' monthlabel ' snow depth']);

% Mean month soil temp by SWE - 20cm
subplot 223;
plot(monthMeans(wasatchTest,7), monthMeans(wasatchTest,4), 'om');
hold on
plot(monthMeans(uintaTest,7), monthMeans(uintaTest,4), 'ob');
x = monthMeans(bothTest, 7);
y = monthMeans(bothTest, 4);
yNoNan = ~isnan(y);
coefficients = polyfit(x(yNoNan), y(yNoNan), 2);
pnFit = polyval(coefficients, (0:1:700));
plot((0:1:700), pnFit,':k');
xlabel('Mean SWE');
ylabel('Mean 20cm soil temp (Celsius)');
title([monthlabel ' 20cm soil temp vs ' monthlabel ' SWE']);

% Mean month soil temp by snow depth - 20cm
subplot 224;
plot(monthMeans(wasatchTest,6), monthMeans(wasatchTest,4), 'om');
hold on
plot(monthMeans(uintaTest,6), monthMeans(uintaTest,4), 'ob');
xlabel('Mean snow depth');
%ylabel('Mean soil temp (Celsius)');
title(['vs ' monthlabel ' snow depth']);

%----------------------------------------------------------------------
% FIG 4 - SoilT vs Air T
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name',[monthlabel ' air vs soil temps - Wasatch and Uinta']);

% Mean month soil temps vs air temps
subplot 221;
plot(monthMeans(wasatchTest, 8), monthMeans(wasatchTest, 3), 'om');
hold on
plot(monthMeans(uintaTest,8), monthMeans(uintaTest, 3), 'ob');
xlabel('Mean AirT');
ylabel('Mean SoilT');
title([monthlabel ' 5cm Ts vs AirT']);

subplot 222;
plot(monthMeans(wasatchTest,8), monthMeans(wasatchTest, 4), 'om');
hold on
plot(monthMeans(uintaTest,8), monthMeans(uintaTest, 4), 'ob');
xlabel('Mean AirT');
ylabel('Mean SoilT');
title([monthlabel ' 20cm Ts vs AirT']);

%----------------------------------------------------------
% FIG 5 - Offsets between AirT and SoilT
offset = (monthMeans(:, 3) - monthMeans(:, 8));

fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name',['Air and soil T offset in ' monthlabel]);

% Mean soil temps vs air temps
subplot 221;
plot(monthMeans(:, 7), offset(:), 'om');
% hold on
% plot(monthMeans(uintaTest,8), monthMeans(uintaTest, 3), 'ob');
xlabel('SWE');
ylabel('AirT-5cmSoilT');
title([monthlabel ' Temperature offset vs SWE']);

subplot 222;
plot(monthMeans(:, 6), offset(:), 'om');
% hold on
% plot(monthMeans(uintaTest,8), monthMeans(uintaTest, 4), 'ob');
xlabel('Snow Depth');
ylabel('AirT-5cmSoilT');
title([monthlabel ' Temperature offset vs Snow Depth']);

% Mean month soil temps vs air temps
subplot 223;
plot(monthMeans(wasatchTest, 7), offset(wasatchTest, 1), 'om');
hold on
plot(monthMeans(uintaTest,7), offset(uintaTest, 1), 'ob');
xlabel('SWE');
ylabel('AirT-5cmSoilT');
legend('Wasatch', 'Uintas');
title([monthlabel ' Temperature offset vs SWE']);

subplot 224;
plot(monthMeans(wasatchTest, 6), offset(wasatchTest, 1), 'om');
hold on
plot(monthMeans(uintaTest, 6), offset(uintaTest, 1), 'ob');
xlabel('Snow Depth');
ylabel('AirT-5cmSoilT');
title([monthlabel ' Temperature offset vs Snow Depth']);

junk = 99;
