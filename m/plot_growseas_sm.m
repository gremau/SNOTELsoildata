function n = plot_growseas_sm()
% plot_growseas_sm.m
%
% Reads the outputs from summarize_wateryear.m and makes a number of plots
% characterizing variability in soil moisture and soil temperature in the
% SNOTEL network.
%
% Feb 20, 2012 - Greg Maurer

clear;          % clear memory
close all;      % clear any figures
fignum=0;       % used to increment figure number for plots

% Add any needed tools
addpath('/home/greg/data/code_resources/m_common/'); 
addpath('/home/greg/data/code_resources/m_common/nanstuff/');
addpath('/home/greg/data/code_resources/m_common/linear/'); 
addpath('/home/greg/data/code_resources/m_common/hline_vline/'); 

% Set data path and file name, read in file
rawdatapath = '../rawdata/soilsensors_hourly/';
processeddatapath = '../processed_data/';

climData = csvread([processeddatapath 'wyear_climatesummary.txt'], 1,0);
vwcData = csvread([processeddatapath 'wyear_soilwatersummary_hourly.txt'], 1,0);
vwcDataN = csvread([processeddatapath ...
    'wyear_soilwatersummary_hourly_smnorm.txt'], 1,0);
soiltempdata = csvread([processeddatapath 'wyear_soiltempsummary_hourly.txt'], 1,0);
soilsiteyears = dlmread([rawdatapath 'filelist.txt'], ',');
soilsites = unique(soilsiteyears(:, 1));

% Parse climate data
site_cl = climData(:, 1);
year_cl = climData(:, 2);
maxswe = climData(:, 3)*25.4;
maxsweday = climData(:, 4);
maxdepth = climData(:, 5);
onsetdoy = climData(:, 6);
meltdoy = climData(:, 7);
snowduration = climData(:, 8);
totaldaysSC = climData(:, 9); % Total days, may vary from duration above
maxcontinSC = climData(:, 10);% Length of longest continuos snowpack
numcontinSC = climData(:, 11);% # of continuous snowcovered periods
accumprecip = climData(:, 12)*25.4;
JASprecip = climData(:, 13);
octSWEmean = climData(:, 19);
octSWEmed = climData(:, 20);
octSWEsd = climData(:, 21);
novSWEmean = climData(:, 22);
novSWEmed = climData(:, 23);
novSWEsd = climData(:, 24);
decSWEmean = climData(:, 25);
decSWEmed = climData(:, 26);
decSWEsd = climData(:, 27);
janSWEmean = climData(:, 28);
janSWEmed = climData(:, 29);
janSWEsd = climData(:, 30);
febSWEmean = climData(:, 31);
febSWEmed = climData(:, 32);
febSWEsd = climData(:, 33);
marSWEmean = climData(:, 34);
marSWEmed = climData(:, 35);
marSWEsd = climData(:, 36);
aprSWEmean = climData(:, 37);
aprSWEmed = climData(:, 38);
aprSWEsd = climData(:, 39);
maySWEmean = climData(:, 40);
maySWEmed = climData(:, 41);
maySWEsd = climData(:, 42);
junSWEmean = climData(:, 43);
junSWEmed = climData(:, 44);
junSWEsd = climData(:, 45);
julSWEmean = climData(:, 46);
julSWEmed = climData(:, 47);
julSWEsd = climData(:, 48);

octTairMean = climData(:, 55);
octTairSd = climData(:, 56);
novTairMean = climData(:, 57);
novTairSd = climData(:, 58);
decTairMean = climData(:, 59);
decTairSd = climData(:, 60);
janTairMean = climData(:, 61);
janTairSd = climData(:, 62);
febTairMean = climData(:, 63);
febTairSd = climData(:, 64);
marTairMean = climData(:, 65);
marTairSd = climData(:, 66);
aprTairMean = climData(:, 67);
aprTairSd = climData(:, 68);
mayTairMean = climData(:, 69);
mayTairSd = climData(:, 70);
junTairMean = climData(:, 71);
junTairSd = climData(:, 72);
julTairMean = climData(:, 73);
julTairSd = climData(:, 74);
augTairMean = climData(:, 75);
augTairSd = climData(:, 76);
sepTairMean = climData(:, 77);
sepTairSd = climData(:, 78);
freemat = climData(:, 79);
scovmat = climData(:, 80);
maat = climData(:, 81);
sdAnnTair = climData(:, 82);
preonsetTair = climData(:, 83);
preonsetTairSd = climData(:, 84);
premeltTair = climData(:, 85);
premeltTairSd = climData(:, 86);
postmeltTair = climData(:, 87);
postmeltTairSd = climData(:, 88);
elev = climData(:, 89);
lat = climData(:, 90);
lon = climData(:, 91);
ltMeanSWE = climData(:, 92);
ltMeanPrecip = climData(:, 93);

% Parse soilwatersummary
site_sw = vwcDataN(:, 1);
year_sw = vwcDataN(:, 2);
octVWC5mean = vwcDataN(:, 3);
octVWC5sd = vwcDataN(:, 4);
octVWC20mean = vwcDataN(:, 5);
octVWC20sd = vwcDataN(:, 6);
octVWC50mean = vwcDataN(:, 7);
octVWC50sd = vwcDataN(:, 8);
novVWC5mean = vwcDataN(:, 9);
novVWC5sd = vwcDataN(:, 10);
novVWC20mean = vwcDataN(:, 11);
novVWC20sd = vwcDataN(:, 12);
novVWC50mean = vwcDataN(:, 13);
novVWC50sd = vwcDataN(:, 14);
decVWC5mean = vwcDataN(:, 15);
decVWC5sd = vwcDataN(:, 16);
decVWC20mean = vwcDataN(:, 17);
decVWC20sd = vwcDataN(:, 18);
decVWC50mean = vwcDataN(:, 19);
decVWC50sd = vwcDataN(:, 20);
janVWC5mean = vwcDataN(:, 21);
janVWC5sd = vwcDataN(:, 22);
janVWC20mean = vwcDataN(:, 23);
janVWC20sd = vwcDataN(:, 24);
janVWC50mean = vwcDataN(:, 25);
janVWC50sd = vwcDataN(:, 26);
febVWC5mean = vwcDataN(:, 27);
febVWC5sd = vwcDataN(:, 28);
febVWC20mean = vwcDataN(:, 29);
febVWC20sd = vwcDataN(:, 30);
febVWC50mean = vwcDataN(:, 31);
febVWC50sd = vwcDataN(:, 32);
marVWC5mean = vwcDataN(:, 33);
marVWC5sd = vwcDataN(:, 34);
marVWC20mean = vwcDataN(:, 35);
marVWC20sd = vwcDataN(:, 36);
marVWC50mean = vwcDataN(:, 37);
marVWC50sd = vwcDataN(:, 38);
aprVWC5mean = vwcDataN(:, 39);
aprVWC5sd = vwcDataN(:, 40);
aprVWC20mean = vwcDataN(:, 41);
aprVWC20sd = vwcDataN(:, 42);
aprVWC50mean = vwcDataN(:, 43);
aprVWC50sd = vwcDataN(:, 44);
mayVWC5mean = vwcDataN(:, 45);
mayVWC5sd = vwcDataN(:, 46);
mayVWC20mean = vwcDataN(:, 47);
mayVWC20sd = vwcDataN(:, 48);
mayVWC50mean = vwcDataN(:, 49);
mayVWC50sd = vwcDataN(:, 50);
junVWC5mean = vwcDataN(:, 51);
junVWC5sd = vwcDataN(:, 52);
junVWC20mean = vwcDataN(:, 53);
junVWC20sd = vwcDataN(:, 54);
junVWC50mean = vwcDataN(:, 55);
junVWC50sd = vwcDataN(:, 56);
julVWC5mean = vwcDataN(:, 57);
julVWC5sd = vwcDataN(:, 58);
julVWC20mean = vwcDataN(:, 59);
julVWC20sd = vwcDataN(:, 60);
julVWC50mean = vwcDataN(:, 61);
julVWC50sd = vwcDataN(:, 62);
augVWC5mean = vwcDataN(:, 63);
augVWC5sd = vwcDataN(:, 64);
augVWC20mean = vwcDataN(:, 65);
augVWC20sd = vwcDataN(:, 66);
augVWC50mean = vwcDataN(:, 67);
augVWC50sd = vwcDataN(:, 68);
sepVWC5mean = vwcDataN(:, 69);
sepVWC5sd = vwcDataN(:, 70);
sepVWC20mean = vwcDataN(:, 71);
sepVWC20sd = vwcDataN(:, 72);
sepVWC50mean = vwcDataN(:, 73);
sepVWC50sd = vwcDataN(:, 74);

% These repeat through sept (end of wy)
ondVWC5mean = vwcData(:, 75);
ondVWC5sd = vwcData(:, 76);
ondVWC20mean = vwcData(:, 77);
ondVWC20sd = vwcData(:, 78);
ondVWC50mean = vwcData(:, 79);
ondVWC50sd = vwcData(:, 81);
jfmVWC5mean = vwcData(:, 81);
jfmVWC5sd = vwcData(:, 82);
jfmVWC20mean = vwcData(:, 83);
jfmVWC20sd = vwcData(:, 84);
jfmVWC50mean = vwcData(:, 85);
jfmVWC50sd = vwcData(:, 86);
amjVWC5mean = vwcData(:, 87);
amjVWC5sd = vwcData(:, 88);
amjVWC20mean = vwcData(:, 89);
amjVWC20sd = vwcData(:, 90);
amjVWC50mean = vwcData(:, 91);
amjVWC50sd = vwcData(:, 92);
jasVWC5mean = vwcData(:, 93);
jasVWC5sd = vwcData(:, 94);
jasVWC20mean = vwcData(:, 95);
jasVWC20sd = vwcData(:, 96);
jasVWC50mean = vwcData(:, 97);
jasVWC50sd = vwcData(:, 98);

% preonsetTair = vwcData(:, 97);
% preonsetVWC5 = vwcData(:, 98);
% preonsetVWC20 = vwcData(:, 99);
% preonsetVWC50 = vwcData(:, 100);
% premeltTair = vwcData(:, 101);
% postmeltTair = vwcData(:, 102);

% Parse NORMALIZED soilwatersummary
amjVWC5meanN = vwcDataN(:, 87);
amjVWC5sdN = vwcDataN(:, 88);
amjVWC20meanN = vwcDataN(:, 89);
amjVWC20sdN = vwcDataN(:, 90);
amjVWC50meanN = vwcDataN(:, 91);
amjVWC50sdN = vwcDataN(:, 92);
jasVWC5meanN = vwcDataN(:, 93);
jasVWC5sdN = vwcDataN(:, 94);
jasVWC20meanN = vwcDataN(:, 95);
jasVWC20sdN = vwcDataN(:, 96);
jasVWC50meanN = vwcDataN(:, 97);
jasVWC50sdN = vwcDataN(:, 98);


% Get a subset of climData that corresponds with available soildata
matchtest = ismember(climData(:, 1:2), soilsiteyears(:, 1:2), 'rows');
matchsets = climData(matchtest, :);

%----------------------------------------------------
% FIG 1 - Plot soil moisture vs max swe snowmelt day
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Wateryear seasonal soil moisture');

subplot (2, 2, 1)
plot(maxswe(matchtest), amjVWC5mean,  '.g', ...
    maxswe(matchtest), amjVWC20mean, '.b', ...
    maxswe(matchtest), amjVWC50mean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Peak SWE (mm)');
ylabel('VWC (%)');
title(' Apr, May, Jun VWC');

subplot (2, 2, 2)
plot(maxswe(matchtest), jasVWC5mean,  '.g', ...
    maxswe(matchtest), jasVWC20mean, '.b', ...
    maxswe(matchtest), jasVWC50mean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Peak SWE (mm)');
ylabel('VWC (%)');
title(' July, Aug, Sept VWC');

subplot (2, 2, 3)
plot(meltdoy(matchtest), amjVWC5mean,  '.g', ...
    meltdoy(matchtest), amjVWC20mean, '.b', ...
    meltdoy(matchtest), amjVWC50mean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Day of wateryear');
ylabel('VWC (%)');
title(' Apr, May, Jun VWC vs snowmelt day');

subplot (2, 2, 4)
plot(meltdoy(matchtest), jasVWC5mean,  '.g', ...
    meltdoy(matchtest), jasVWC20mean, '.b', ...
    meltdoy(matchtest), jasVWC50mean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Day of wateryear (mm)');
ylabel('VWC (%)');
title(' July, Aug, Sept VWC vs snowmeltday');

%----------------------------------------------------
% FIG 2 - Plot NORMALIZED soil moisture vs max swe  and snowmelt day
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Wateryear NORMALIZED seasonal soil moisture');

subplot (2, 2, 1)
plot(maxswe(matchtest), amjVWC5meanN,  '.g', ...
    maxswe(matchtest), amjVWC20meanN, '.b', ...
    maxswe(matchtest), amjVWC50meanN, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Peak SWE (mm)');
ylabel('VWC (%)');
title(' Apr, May, Jun VWC');

subplot (2, 2, 2)
plot(maxswe(matchtest), jasVWC5meanN,  '.g', ...
    maxswe(matchtest), jasVWC20meanN, '.b', ...
    maxswe(matchtest), jasVWC50meanN, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Peak SWE (mm)');
ylabel('VWC (%)');
title(' July, Aug, Sept VWC');

subplot (2, 2, 3)
plot(meltdoy(matchtest), amjVWC5meanN,  '.g', ...
    meltdoy(matchtest), amjVWC20meanN, '.b', ...
    meltdoy(matchtest), amjVWC50meanN, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Day of wateryear');
ylabel('VWC (%)');
title(' Apr, May, Jun VWC vs snowmelt day');

subplot (2, 2, 4)
plot(meltdoy(matchtest), jasVWC5meanN,  '.g', ...
    meltdoy(matchtest), jasVWC20meanN, '.b', ...
    meltdoy(matchtest), jasVWC50meanN, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Day of wateryear (mm)');
ylabel('VWC (%)');
title(' July, Aug, Sept VWC vs snowmeltday');


%----------------------------------------------------
% FIG 3 - Plot soil moisture vs snowmelt day and peak SWE
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Jul, Aug, Sep VWC vs snowpack');

% Left side - plot vs Snowmelt day
% Set binning parameters
topEdge = 305; % define limits
botEdge = 145; % define limits
numBins = 15; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1))+5;

x = meltdoy(matchtest); %split into x and y
y = jasVWC5mean;
y2 = jasVWC5sd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(321); % 5cm  
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([140 310]); ylim([0 45]);
title('Summer VWC vs snowmelt day'); 
ylabel('5cm VWC');

y = jasVWC20mean;
y2 = jasVWC20sd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(323); % 20cm  
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([140 310]); ylim([0 45]);
ylabel('20cm VWC');

% 50cm
y = jasVWC50mean;
y2 = jasVWC50sd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(325); % 50cm  
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([140 310]); ylim([0 45]);
xlabel('Snowmelt day'); ylabel('50cm VWC');

% Right side - plot vs max SWE
topEdge = 1600; % define limits
botEdge = 0; % define limits
numBins = 15; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1))+53;

x = maxswe(matchtest); %split into x and y
y = jasVWC5mean;
y2 = jasVWC5sd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(322);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([0 1650]); ylim([0 45]);
title('Summer VWC vs peak SWE');
ylabel('5cm VWC');

y = jasVWC20mean;
y2 = jasVWC20sd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(324);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([0 1650]); ylim([0 45]);
ylabel('20cm VWC');

y = jasVWC50mean;
y2 = jasVWC50sd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(326);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([0 1650]); ylim([0 45]);
xlabel('Peak SWE (mm)');

%----------------------------------------------------
% FIG 4 - Plot NORMALIZED soil moisture vs snowmelt day and peak SWE
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Normalized Jul, Aug, Sep VWC vs snowpack');

% Left side - plot vs Snowmelt day
% Set binning parameters
topEdge = 305; % define limits
botEdge = 145; % define limits
numBins = 15; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1))+5;

x = meltdoy(matchtest); %split into x and y
y = jasVWC5meanN;
y2 = jasVWC5sdN;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(321); % 5cm  
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([140 310]); ylim([0 1]);
title('Summer VWC vs snowmelt day'); 
ylabel('5cm VWC-norm');

y = jasVWC20meanN;
y2 = jasVWC20sdN;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(323); % 20cm  
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([140 310]); ylim([0 1]);
ylabel('20cm VWC-norm');

% 50cm
y = jasVWC50meanN;
y2 = jasVWC50sdN;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(325); % 50cm  
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([140 310]); ylim([0 1]);
xlabel('Snowmelt day'); ylabel('50cm VWC-norm');

% Right side - plot vs max SWE
topEdge = 1600; % define limits
botEdge = 0; % define limits
numBins = 15; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1))+53;

x = maxswe(matchtest); %split into x and y
y = jasVWC5meanN;
y2 = jasVWC5sdN;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(322);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([0 1650]); ylim([0 1]);
title('Summer VWC vs peak SWE'); 
ylabel('5cm VWC');

y = jasVWC20meanN;
y2 = jasVWC20sdN;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(324);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([0 1650]); ylim([0 1]);
ylabel('20cm VWC');

y = jasVWC50meanN;
y2 = jasVWC50sdN;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(326);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([0 1650]); ylim([0 1]);
xlabel('Peak SWE (mm)');

%----------------------------------------------------
% FIG 5 - Plot NORMALIZED soil moisture vs snowmelt day and peak SWE
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Normalized Jul, Aug, Sep VWC vs snowpack');

polyorder = 1;

x = meltdoy(matchtest); %split into x and y
y = jasVWC5meanN;
subplot(321); % 5cm  
errorbar(x, y, jasVWC5sdN, 'o', 'Color', [0.8 0.8 0.8],...
    'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', 'w');
hold on;
[~, rsq, xfit, yfit] = fitline(x, y, polyorder, [100, 400]);
plot(xfit, yfit, '-k', 'LineWidth', 2);
text(150, 0.75, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlim([140 310]); ylim([0 1]);
title('Summer VWC vs snowmelt day'); 
ylabel('5cm VWC-norm');

y = jasVWC20meanN;
subplot(323); % 20cm  
errorbar(x, y, jasVWC20sdN, 'o', 'Color', [0.8 0.8 0.8],...
    'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', 'w');
hold on;
[~, rsq, xfit, yfit] = fitline(x, y, polyorder, [100, 400]);
plot(xfit, yfit, '-k', 'LineWidth', 2);
text(150, 0.75, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlim([140 310]); ylim([0 1]);
ylabel('20cm VWC-norm');

y = jasVWC50meanN;
subplot(325); % 50cm  
errorbar(x, y, jasVWC50sdN, 'o', 'Color', [0.8 0.8 0.8],...
    'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', 'w');
hold on
[~, rsq, xfit, yfit] = fitline(x, y, polyorder, [100, 400]);
plot(xfit, yfit, '-k', 'LineWidth', 2);
text(150, 0.75, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlim([140 310]); ylim([0 1]);
ylabel('50cm VWC-norm');

% Right side - plot vs max SWE

x = maxswe(matchtest); %split into x and y
y = jasVWC5meanN;
subplot(322); % 5cm  
errorbar(x, y, jasVWC5sdN, 'o', 'Color', [0.8 0.8 0.8],...
    'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', 'w');
hold on;
[~, rsq, xfit, yfit] = fitline(x, y, polyorder, [0, 1700]);
plot(xfit, yfit, '-k', 'LineWidth', 2);
text(150, 0.75, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlim([0 1650]); ylim([0 1]);
title('Summer VWC vs peak SWE'); 

y = jasVWC20meanN;
subplot(324); % 20cm  
errorbar(x, y, jasVWC20sdN, 'ok', 'Color', [0.8 0.8 0.8],...
    'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', 'w');
hold on;
[~, rsq, xfit, yfit] = fitline(x, y, polyorder, [0, 1700]);
plot(xfit, yfit, '-k', 'LineWidth', 2);
text(150, 0.75, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlim([0 1650]); ylim([0 1]);

y = jasVWC50meanN;
subplot(326); % 50cm  
errorbar(x, y, jasVWC50sdN, 'o', 'Color', [0.8 0.8 0.8],...
    'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', 'w');
hold on
[~, rsq, xfit, yfit] = fitline(x, y, polyorder, [0, 1700]);
plot(xfit, yfit, '-k', 'LineWidth', 2);
text(150, 0.75, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlim([0 1650]); ylim([0 1]);
xlabel('Peak SWE (mm)');

%----------------------------------------------------
% FIG 6 - Plot August soil moisture vs snowmelt date
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Aug soil moisture vs snomelt date');

% Left side - plot vs Snowmelt day
% Set binning parameters
topEdge = 305; % define limits
botEdge = 175; % define limits
numBins = 12; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1))+5;

% 5cm
x = meltdoy(matchtest); %split into x and y
y = augVWC5mean;
y2 = augVWC5sd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(221);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([170 310]); %ylim([0 1]);
ticks = [175;240;300]
set(gca,'XTick', ticks, 'XTickLabel', {'Apr';'Jun'; 'Aug'});
xlabel('Snowmelt dowy');
title('5cm VWC');

y = augVWC20mean;
y2 = augVWC20sd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(222);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([170 310]);%ylim([0 1]);
ticks = [180;240;300]
set(gca,'XTick', ticks, 'XTickLabel', {'Apr';'Jun'; 'Aug'});
xlabel('Snowmelt dowy');
title('20cm VWC');

y = augVWC50mean;
y2 = augVWC50sd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(223);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([170 310]);%ylim([0 1]);
ticks = [180;240;300]
set(gca,'XTick', ticks, 'XTickLabel', {'Apr';'Jun'; 'Aug'});
xlabel('Snowmelt dowy');
ylabel('Aug VWC (normalized)')
title('50cm VWC');

%--------------------------------------------------------------
% FIG 7 - Regress snowpack vs growing season vwc for 3 sites
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Regress vwc vs snowpack for 3 sites');
test = site_cl==828;
site1meltd = meltdoy(matchtest & test);
site1pks = maxswe(matchtest & test);
site1sm20 = jasVWC20meanN(site_sw==828);
test = site_cl==330;
site2meltd = meltdoy(matchtest & test);
site2pks = maxswe(matchtest & test);
site2sm20 = jasVWC20meanN(site_sw==330);
test = site_cl==582;
site3meltd = meltdoy(matchtest & test);
site3pks = maxswe(matchtest & test);
site3sm20 = jasVWC20meanN(site_sw==582);

subplot (3, 2, 1)
plot(site1meltd, site1sm20, '.b');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site1meltd, site1sm20, 1, xrange);
plot(xfit, yfit,':k');
text(240, .8, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('Snowmelt Day'); ylabel('20cm VWC');
title('Trial Lake');

subplot (3, 2, 3)
plot(site2meltd, site2sm20, '.b');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site2meltd, site2sm20, 1, xrange);
plot(xfit, yfit,':k');
text(210, .2, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('Snowmelt Day'); ylabel('20cm VWC');
title('Beaver Divide');

subplot (3, 2, 5)
plot(site3meltd, site3sm20, '.b');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site3meltd, site3sm20, 1, xrange);
plot(xfit, yfit,':k');
text(210, .45, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('Snowmelt Day'); ylabel('20cm VWC');
title('Little Bear');

subplot (3, 2, 2)
plot(site1pks, site1sm20, '.b')
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site1pks, site1sm20, 1, xrange);
plot(xfit, yfit,':k');
text(200, .7, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('Peak SWE'); ylabel('20cm VWC');
title('Trial Lake');

subplot (3, 2, 4)
plot(site2pks, site2sm20, '.b')
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site2pks, site2sm20, 1, xrange);
plot(xfit, yfit,':k');
text(350, .2, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('Peak SWE'); ylabel('20cm VWC');
title('Beaver Divide');

subplot (3, 2, 6)
plot(site3pks, site3sm20, '.b')
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site3pks, site3sm20, 1, xrange);
plot(xfit, yfit,':k');
text(250, .45, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('Peak SWE'); ylabel('20cm VWC');
title('Little Bear');

% -----------------------------------------------------------
% FIG 4 - Month soil water content vs elevation
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name', 'Mean monthly 20cm VWC vs elevation - All SNOTEL/years');

% Set some plotting parameters
plotorder = 1:3;
months = ['oct';'nov';'dec';'jan';'feb';'mar';'apr';'may';...
    'jun';'jul';'aug';'sep'] ;
polyorder = 1;

for i=plotorder
    subplot(4, 3, i);
    x = elev(matchtest);
    eval(['y1 = ' months(i,:) 'VWC20mean;']);
    %eval(['y2 = ' months(i,:) 'Ts5mean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    [~, rsq, xfit, yfit] = fitline(x, y1, polyorder, [800, 3600]);
    plot(xfit, yfit, '-k', 'LineWidth', 2);
    text(0.1, 0.9, ['r^2 = ' num2str(rsq, 2)], 'units', 'normalized'); % r^2 values
    %plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    title(months(i,:));
    xlim([800, 3600]); ylim([0, 1]);
    hline(0, '--k');
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('^oC')
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
    
    subplot(4,3, i+3)
    eval(['y1 = ' months(i+3,:) 'VWC20mean;']);
    %eval(['y2 = ' months(i+3,:) 'Ts20mean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    [~, rsq, xfit, yfit] = fitline(x, y1, polyorder, [800, 3600]);
    plot(xfit, yfit, '-k', 'LineWidth', 2);
    text(0.1, 0.9, ['r^2 = ' num2str(rsq, 2)], 'units', 'normalized'); % r^2 values
%     plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    title(months(i+3,:));
    xlim([800, 3600]); ylim([0, 1]);
    hline(0, '--k');
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('^oC')
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
    
    subplot(4,3, i+6);
    eval(['y1 = ' months(i+6,:) 'VWC20mean;']);
    %eval(['y2 = ' months(i+6,:) 'Ts50mean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    [~, rsq, xfit, yfit] = fitline(x, y1, polyorder, [800, 3600]);
    plot(xfit, yfit, '-k', 'LineWidth', 2);
    text(0.1, 0.9, ['r^2 = ' num2str(rsq, 2)], 'units', 'normalized'); % r^2 values
%     plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    title(months(i+6,:));
    xlim([800, 3600]); ylim([0, 1]);
    hline(0, '--k');
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('^oC')
    elseif i>1
        set(gca, 'YTickLabel', '');
    end

    subplot(4,3, i+9);
    eval(['y1 = ' months(i+9,:) 'VWC20mean;']);
    %eval(['y2 = ' months(i+9,:) 'Ts50mean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    [~, rsq, xfit, yfit] = fitline(x, y1, polyorder, [800, 3600]);
    plot(xfit, yfit, '-k', 'LineWidth', 2);
    text(0.1, 0.9, ['r^2 = ' num2str(rsq, 2)], 'units', 'normalized'); % r^2 values
%     plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    title(months(i+9,:));
    xlim([800, 3600]); ylim([0, 1]);
    hline(0, '--k');
    if i==1
        ylabel('^oC')
    elseif i==2
        xlabel('Elevation (m)');
        set(gca, 'YTickLabel', '');
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
end


