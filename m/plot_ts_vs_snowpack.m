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
addpath('~/data/code_resources/m_common/nanstuff/');

% Ask user for month number
monthsel = str2double(input('Which month (1-12)?: ', 's'));

% Set data path and file name, read in file
rawdatapath = '../rawdata/';
processeddatapath = '../processed_data/';

% Load list of sites with data in the daily data directory
dailysites = sortrows(csvread('../rawdata/allsensors_daily/filelist.txt'));
soilsites = sortrows(csvread('../rawdata/soilsensors_hourly/filelist.txt'));

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

% LOAD the data (can switch between daily/hourly data here
climData = csvread([processeddatapath 'wyear_climatesummary.txt']);
% soilsites = unique(filelistSoil(:, 1));
% Soil temp data
tsData = csvread([processeddatapath 'wyear_soiltempsummary_hourly.txt']);
% tsData = csvread([processeddatapath 'wyear_soiltempsummary_daily.txt']);

% Get a subset of climData that corresponds with available soildata
[matchtest, idx] = ismember(climData(:, 1:2), tsData(:, 1:2), 'rows');
soilClim = climData(matchtest, :);

% Now assign variables
octSWEmean = climData(:, 14);
octSWEmed = climData(:, 15);
octSWEsd = climData(:, 16);
novSWEmean = climData(:, 17);
novSWEmed = climData(:, 18);
novSWEsd = climData(:, 19);
decSWEmean = climData(:, 20);
decSWEmed = climData(:, 21);
decSWEsd = climData(:, 22);
janSWEmean = climData(:, 23);
janSWEmed = climData(:, 24);
janSWEsd = climData(:, 25);
febSWEmean = climData(:, 26);
febSWEmed = climData(:, 27);
febSWEsd = climData(:, 28);
marSWEmean = climData(:, 29);
marSWEmed = climData(:, 30);
marSWEsd = climData(:, 31);
aprSWEmean = climData(:, 32);
aprSWEmed = climData(:, 33);
aprSWEsd = climData(:, 34);
maySWEmean = climData(:, 35);
maySWEmed = climData(:, 36);
maySWEsd = climData(:, 37);
junSWEmean = climData(:, 38);
junSWEmed = climData(:, 39);
junSWEsd = climData(:, 40);
julSWEmean = climData(:, 41);
julSWEmed = climData(:, 42);
julSWEsd = climData(:, 43);

oct5cmSTmean = tsData(:, 3);
oct5cmSTsd = tsData(:, 4);
oct20cmSTmean = tsData(:, 5);
oct20cmSTsd = tsData(:, 6);
oct50cmSTmean = tsData(:, 7);
oct50cmSTsd = tsData(:, 8);
dec5cmSTmean = tsData(:, 15);
dec5cmSTsd = tsData(:, 16);
dec20cmSTmean = tsData(:, 17);
dec20cmSTsd = tsData(:, 18);
dec50cmSTmean = tsData(:, 19);
dec50cmSTsd = tsData(:, 20);
% These repeat through sept (end of wy)


% sites = unique(soilsites(:, 1));
% monthLabels = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sept' 'Oct'...
%     'Nov' 'Dec'};
% monthlabel = monthLabels{monthsel};
% d_monthMeans = [];
% h_monthMeans = [];
% % Load data and parse out month data
% for i = 1:length(sites);
%     dData = loadsnotel(sites(i), 'daily', 'exclude');
%     hData = loadsnotel(sites(i), 'hourly', 'exclude');
%     %Create datevector for datafile
%     dDatevec = datevec(dData{2}, 'yyyy-mm-dd');
%     hDatevec = datevec(strcat(hData{2}, hData{3}), 'yyyy-mm-ddHH:MM');
%     % Columns are site, year, sndepth, swe, airT
%     dailyData = [double(dData{1}) dDatevec(:,1) dData{10}*25.4 ...
%         dData{4}*25.4 dData{9}];
%     % Columns are site, year ts5, ts20, ts50
%     hourlyData = [double(hData{1}) hDatevec(:,1) hData{7} hData{8} hData{9}];
%     % Get monthly data
%     d_monthTest = dDatevec(:,2)==monthsel;
%     h_monthTest = hDatevec(:,2)==monthsel;
%     d_monthData = dailyData(d_monthTest, :);
%     h_monthData = hourlyData(h_monthTest, :);
%     monthYears = unique(h_monthData(:,2));
%     % Reduce to yearly averages
%     for j = 1:length(monthYears)
%         yearTest1 = d_monthData(:, 2) == monthYears(j);
%         yearTest2 = h_monthData(:, 2) == monthYears(j);
%         d_monthMeans = [d_monthMeans; nanmean(d_monthData(yearTest1, :), 1)];
%         h_monthMeans = [h_monthMeans; nanmean(h_monthData(yearTest2, :), 1)];
%     end
% end
% 
% % Verify that the site/year rows in d_ and h_monthMeans are the same
% dhlogical = ismember(d_monthMeans(:,1:2), h_monthMeans(:,1:2), 'rows');
% sum(dhlogical)==length(h_monthMeans)
% 
% % site, year ts5, ts20, ts50, sndepth, swe, airT
% monthMeans = [h_monthMeans, d_monthMeans(:, 3:5)];



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
