% plot_ts_snowpack.m
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
addpath('~/data/code_resources/m_common/');
addpath('~/data/code_resources/m_common/nanstuff/');

% Ask user for month number
%monthsel = str2double(input('Which month (1-12)?: ', 's'));

% Set processed data path
processeddatapath = '../processed_data/';

% Load lists of sites with data in the daily/hourly data directory
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
clear wasatchTest uintaTest wasatchUintaCell;

% LOAD the data (can switch between daily/hourly data here)
climData = csvread([processeddatapath 'wyear_climatesummary.txt']);
tsData = csvread([processeddatapath 'wyear_soiltempsummary_hourly.txt']);
% tsData = csvread([processeddatapath 'wyear_soiltempsummary_daily.txt']);

% Get a subset of climData that corresponds with available soildata
[matchtest, idx] = ismember(climData(:, 1:2), tsData(:, 1:2), 'rows');
soilClim = climData(matchtest, :);
% matchtest2 = ismember(tsData(:, 1:2), soilClim(:, 1:2), 'rows');

% Now assign variables
octSWEmean = climData(:, 14)*25.4;
octSWEmed = climData(:, 15)*25.4;
octSWEsd = climData(:, 16)*25.4;
novSWEmean = soilClim(:, 17)*25.4;
novSWEmed = soilClim(:, 18)*25.4;
novSWEsd = soilClim(:, 19)*25.4;
decSWEmean = soilClim(:, 20)*25.4;
decSWEmed = soilClim(:, 21)*25.4;
decSWEsd = soilClim(:, 22)*25.4;
janSWEmean = soilClim(:, 23)*25.4;
janSWEmed = soilClim(:, 24)*25.4;
janSWEsd = soilClim(:, 25)*25.4;
febSWEmean = soilClim(:, 26)*25.4;
febSWEmed = soilClim(:, 27)*25.4;
febSWEsd = soilClim(:, 28)*25.4;
marSWEmean = soilClim(:, 29)*25.4;
marSWEmed = soilClim(:, 30)*25.4;
marSWEsd = soilClim(:, 31)*25.4;
aprSWEmean = soilClim(:, 32)*25.4;
aprSWEmed = soilClim(:, 33)*25.4;
aprSWEsd = soilClim(:, 34)*25.4;
maySWEmean = soilClim(:, 35)*25.4;
maySWEmed = soilClim(:, 36)*25.4;
maySWEsd = soilClim(:, 37)*25.4;
junSWEmean = soilClim(:, 38)*25.4;
junSWEmed = soilClim(:, 39)*25.4;
junSWEsd = soilClim(:, 40)*25.4;
julSWEmean = soilClim(:, 41)*25.4;
julSWEmed = soilClim(:, 42)*25.4;
julSWEsd = soilClim(:, 43)*25.4;

oct5cmSTmean = tsData(:, 3);
oct5cmSTsd = tsData(:, 4);
oct20cmSTmean = tsData(:, 5);
oct20cmSTsd = tsData(:, 6);
oct50cmSTmean = tsData(:, 7);
oct50cmSTsd = tsData(:, 8);
nov5cmSTmean = tsData(:, 9);
nov5cmSTsd = tsData(:, 10);
nov20cmSTmean = tsData(:, 11);
nov20cmSTsd = tsData(:, 12);
nov50cmSTmean = tsData(:, 13);
nov50cmSTsd = tsData(:, 14);
dec5cmSTmean = tsData(:, 15);
dec5cmSTsd = tsData(:, 16);
dec20cmSTmean = tsData(:, 17);
dec20cmSTsd = tsData(:, 18);
dec50cmSTmean = tsData(:, 19);
dec50cmSTsd = tsData(:, 20);
jan5cmSTmean = tsData(:, 21);
jan5cmSTsd = tsData(:, 22);
jan20cmSTmean = tsData(:, 23);
jan20cmSTsd = tsData(:, 24);
jan50cmSTmean = tsData(:, 25);
jan50cmSTsd = tsData(:, 26);
feb5cmSTmean = tsData(:, 27);
feb5cmSTsd = tsData(:, 28);
feb20cmSTmean = tsData(:, 29);
feb20cmSTsd = tsData(:, 30);
feb50cmSTmean = tsData(:, 31);
feb50cmSTsd = tsData(:, 32);
mar5cmSTmean = tsData(:, 33);
mar5cmSTsd = tsData(:, 34);
mar20cmSTmean = tsData(:, 35);
mar20cmSTsd = tsData(:, 36);
mar50cmSTmean = tsData(:, 37);
mar50cmSTsd = tsData(:, 38);
apr5cmSTmean = tsData(:, 39);
apr5cmSTsd = tsData(:, 40);
apr20cmSTmean = tsData(:, 41);
apr20cmSTsd = tsData(:, 42);
apr50cmSTmean = tsData(:, 43);
apr50cmSTsd = tsData(:, 44);
may5cmSTmean = tsData(:, 45);
may5cmSTsd = tsData(:, 46);
may20cmSTmean = tsData(:, 47);
may20cmSTsd = tsData(:, 48);
may50cmSTmean = tsData(:, 49);
may50cmSTsd = tsData(:, 50);
jun5cmSTmean = tsData(:, 51);
jun5cmSTsd = tsData(:, 52);
jun20cmSTmean = tsData(:, 53);
jun20cmSTsd = tsData(:, 54);
jun50cmSTmean = tsData(:, 55);
jun50cmSTsd = tsData(:, 56);
% These repeat through sept (end of wy)


% PLOTS
%----------------------------------------------------------------------
% FIG 1 - Month soil temps vs mean SWE
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name', 'Mean monthly temp vs Mean SWE - all sites/wateryears');

% Set some plotting parameters
plotorder = 1:5;
months = ['dec';'jan';'feb';'mar';'apr'] ;
polyorder = 3;

for i=plotorder
    subplot(3, 5, i);
    eval(['x = ' months(i,:) 'SWEmean;']);
    eval(['y = ' months(i,:) '5cmSTmean;']);
    plot(x, y, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    [~, rsq, xfit, yfit] = fitline(x, y, polyorder, [0, 1500]);
    plot(xfit, yfit, '-k', 'LineWidth', 2);
    text(700, 8, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
    xlim([0, 1500]); ylim([-10, 10]);
    
    subplot(3, 5, i+5)
    eval(['x = ' months(i,:) 'SWEmean;']);
    eval(['y = ' months(i,:) '20cmSTmean;']);
    plot(x, y, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    [~, rsq, xfit, yfit] = fitline(x, y, polyorder, [0, 1500]);
    plot(xfit, yfit, '-k', 'LineWidth', 2);
    text(700, 8, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
    xlim([0, 1500]); ylim([-10, 10]);
    
    subplot(3, 5, i+10);
    eval(['x = ' months(i,:) 'SWEmean;']);
    eval(['y = ' months(i,:) '50cmSTmean;']);
    plot(x, y, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    [~, rsq, xfit, yfit] = fitline(x, y, polyorder, [0, 1500]);
    plot(xfit, yfit, '-k', 'LineWidth', 2);
    text(700, 8, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
    xlim([0, 1500]); ylim([-10, 10]);
end


%----------------------------------------------------------------------
% FIG 2 - Month soil temps vs median SWE
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name', 'Mean monthly temp vs Median SWE - all sites/wateryears');

% plotorder, months, and polyorder are set in previous figure's code

for i=plotorder
    subplot(3, 5, i)
    eval(['x = ' months(i,:) 'SWEmed;']);
    eval(['y = ' months(i,:) '5cmSTmean;']);
    plot(x, y, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    [~, rsq, xfit, yfit] = fitline(x, y, polyorder, [0, 1500]);
    plot(xfit, yfit, '-k', 'LineWidth', 2);
    text(700, 8, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
    xlim([0, 1500]); ylim([-10, 10]);
    
    subplot(3, 5, i+5)
    eval(['x = ' months(i,:) 'SWEmed;']);
    eval(['y = ' months(i,:) '20cmSTmean;']);
    plot(x, y, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    [~, rsq, xfit, yfit] = fitline(x, y, polyorder, [0, 1500]);
    plot(xfit, yfit, '-k', 'LineWidth', 2);
    text(700, 8, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
    xlim([0, 1500]); ylim([-10, 10]);
    
    subplot(3, 5, i+10)
    eval(['x = ' months(i,:) 'SWEmed;']);
    eval(['y = ' months(i,:) '50cmSTmean;']);
    plot(x, y, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    [~, rsq, xfit, yfit] = fitline(x, y, polyorder, [0, 1500]);
    plot(xfit, yfit, '-k', 'LineWidth', 2);
    text(700, 8, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
    xlim([0, 1500]); ylim([-10, 10]);
end

%--------------------------------------------------------
% FIG 2 - Same as above, but tweaked for AGU 2011 poster
fignum = fignum+1;    
h = figure(fignum);

subplot 121;
plot(janSWEmean, jan5cmSTmean, 'ob');
xlabel('Mean SWE (mm)');
ylabel('Soil T (^oC)');
title('December soil temperatures');
legend('5cm one-month mean');

% Mean month soil temp by SWE - 20cm
subplot 122;
plot(janSWEmean, jan20cmSTmean, 'ok');
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
