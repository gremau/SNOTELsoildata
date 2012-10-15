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

% access to nanmean, etc
addpath('/home/greg/data/code_resources/m_common/nanstuff/');

% Set data path and file name, read in file
rawdatapath = '../rawdata/soilsensors_hourly/';
processeddatapath = '../processed_data/';

climData = csvread([processeddatapath 'wyear_climatesummary.txt']);
vwcData = csvread([processeddatapath 'wyear_soilwatersummary_hourly.txt']);
vwcDataN = csvread([processeddatapath ...
    'wyear_soilwatersummary_hourly_smnorm.txt']);
soiltempdata = csvread([processeddatapath 'wyear_soiltempsummary_hourly.txt']);
soilsiteyears = dlmread([rawdatapath 'filelist.txt'], ',');
soilsites = unique(soilsiteyears(:, 1));

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

octAirTmean = climData(:, 50);
octAirTsd = climData(:, 51);
novAirTmean = climData(:, 52);
novAirTsd = climData(:, 53);
decAirTmean = climData(:, 54);
decAirTsd = climData(:, 55);
janAirTmean = climData(:, 56);
janAirTsd = climData(:, 57);
febAirTmean = climData(:, 58);
febAirTsd = climData(:, 59);
marAirTmean = climData(:, 60);
marAirTsd = climData(:, 61);
aprAirTmean = climData(:, 62);
aprAirTsd = climData(:, 63);
mayAirTmean = climData(:, 64);
mayAirTsd = climData(:, 65);
junAirTmean = climData(:, 66);
junAirTsd = climData(:, 67);
julAirTmean = climData(:, 68);
julAirTsd = climData(:, 69);
augAirTmean = climData(:, 70);
augAirTsd = climData(:, 71);
sepAirTmean = climData(:, 72);
sepAirTsd = climData(:, 73);
maat = climData(:, 74);
sdAnnAirT = climData(:, 75);

preonsetAirT = climData(:, 76);
preonsetAirTsd = climData(:, 77);
premeltAirT = climData(:, 78);
premeltAirTsd = climData(:, 79);
postmeltAirT = climData(:, 80);
postmeltAirTsd = climData(:, 81);
elev = climData(:, 82);
lat = climData(:, 83);
lon = climData(:, 84);
ltMeanSWE = climData(:, 85);
ltMeanPrecip = climData(:, 86);

% Parse soilwatersummary
site_sw = vwcData(:, 1);
year_sw = vwcData(:, 2);
oct5cmSMmean = vwcData(:, 3);
oct5cmSMsd = vwcData(:, 4);
oct20cmSMmean = vwcData(:, 5);
oct20cmSMsd = vwcData(:, 6);
oct50cmSMmean = vwcData(:, 7);
oct50cmSMsd = vwcData(:, 8);

feb5cmSMmean = vwcData(:, 27);
feb5cmSMsd = vwcData(:, 28);
feb20cmSMmean = vwcData(:, 29);
feb20cmSMsd = vwcData(:, 30);
feb50cmSMmean = vwcData(:, 31);
feb50cmSMsd = vwcData(:, 32);

may5cmSMmean = vwcData(:, 45);
may5cmSMsd = vwcData(:, 46);
may20cmSMmean = vwcData(:, 47);
may20cmSMsd = vwcData(:, 48);
may50cmSMmean = vwcData(:, 49);
may50cmSMsd = vwcData(:, 50);

aug5cmSMmean = vwcData(:, 63);
aug5cmSMsd = vwcData(:, 64);
aug20cmSMmean = vwcData(:, 65);
aug20cmSMsd = vwcData(:, 66);
aug50cmSMmean = vwcData(:, 67);
aug50cmSMsd = vwcData(:, 68);

% These repeat through sept (end of wy)
ond5cmSMmean = vwcData(:, 75);
ond5cmSMsd = vwcData(:, 76);
ond20cmSMmean = vwcData(:, 77);
ond20cmSMsd = vwcData(:, 78);
ond50cmSMmean = vwcData(:, 79);
ond50cmSMsd = vwcData(:, 81);
jfm5cmSMmean = vwcData(:, 81);
jfm5cmSMsd = vwcData(:, 82);
jfm20cmSMmean = vwcData(:, 83);
jfm20cmSMsd = vwcData(:, 84);
jfm50cmSMmean = vwcData(:, 85);
jfm50cmSMsd = vwcData(:, 86);
amj5cmSMmean = vwcData(:, 87);
amj5cmSMsd = vwcData(:, 88);
amj20cmSMmean = vwcData(:, 89);
amj20cmSMsd = vwcData(:, 90);
amj50cmSMmean = vwcData(:, 91);
amj50cmSMsd = vwcData(:, 92);
jas5cmSMmean = vwcData(:, 93);
jas5cmSMsd = vwcData(:, 94);
jas20cmSMmean = vwcData(:, 95);
jas20cmSMsd = vwcData(:, 96);
jas50cmSMmean = vwcData(:, 97);
jas50cmSMsd = vwcData(:, 98);

% preonsetAirT = vwcData(:, 97);
% preonset5cmSM = vwcData(:, 98);
% preonset20cmSM = vwcData(:, 99);
% preonset50cmSM = vwcData(:, 100);
% premeltAirT = vwcData(:, 101);
% postmeltAirT = vwcData(:, 102);


% Parse NORMALIZED soilwatersummary
amj5cmSMmeanN = vwcDataN(:, 87);
amj5cmSMsdN = vwcDataN(:, 88);
amj20cmSMmeanN = vwcDataN(:, 89);
amj20cmSMsdN = vwcDataN(:, 90);
amj50cmSMmeanN = vwcDataN(:, 91);
amj50cmSMsdN = vwcDataN(:, 92);
jas5cmSMmeanN = vwcDataN(:, 93);
jas5cmSMsdN = vwcDataN(:, 94);
jas20cmSMmeanN = vwcDataN(:, 95);
jas20cmSMsdN = vwcDataN(:, 96);
jas50cmSMmeanN = vwcDataN(:, 97);
jas50cmSMsdN = vwcDataN(:, 98);


% Get a subset of climData that corresponds with available soildata
matchtest = ismember(climData(:, 1:2), soilsiteyears(:, 1:2), 'rows');
matchsets = climData(matchtest, :);


%----------------------------------------------------
% FIG 1 - Plot soil moisture vs max swe snowmelt day
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Wateryear seasonal soil moisture');

subplot (2, 2, 1)
plot(maxswe(matchtest), amj5cmSMmean,  '.g', ...
    maxswe(matchtest), amj20cmSMmean, '.b', ...
    maxswe(matchtest), amj50cmSMmean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Peak SWE (mm)');
ylabel('VWC (%)');
title(' Apr, May, Jun VWC');

subplot (2, 2, 2)
plot(maxswe(matchtest), jas5cmSMmean,  '.g', ...
    maxswe(matchtest), jas20cmSMmean, '.b', ...
    maxswe(matchtest), jas50cmSMmean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Peak SWE (mm)');
ylabel('VWC (%)');
title(' July, Aug, Sept VWC');

subplot (2, 2, 3)
plot(meltdoy(matchtest), amj5cmSMmean,  '.g', ...
    meltdoy(matchtest), amj20cmSMmean, '.b', ...
    meltdoy(matchtest), amj50cmSMmean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Day of wateryear');
ylabel('VWC (%)');
title(' Apr, May, Jun VWC vs snowmelt day');

subplot (2, 2, 4)
plot(meltdoy(matchtest), jas5cmSMmean,  '.g', ...
    meltdoy(matchtest), jas20cmSMmean, '.b', ...
    meltdoy(matchtest), jas50cmSMmean, '.k');
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
plot(maxswe(matchtest), amj5cmSMmeanN,  '.g', ...
    maxswe(matchtest), amj20cmSMmeanN, '.b', ...
    maxswe(matchtest), amj50cmSMmeanN, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Peak SWE (mm)');
ylabel('VWC (%)');
title(' Apr, May, Jun VWC');

subplot (2, 2, 2)
plot(maxswe(matchtest), jas5cmSMmeanN,  '.g', ...
    maxswe(matchtest), jas20cmSMmeanN, '.b', ...
    maxswe(matchtest), jas50cmSMmeanN, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Peak SWE (mm)');
ylabel('VWC (%)');
title(' July, Aug, Sept VWC');

subplot (2, 2, 3)
plot(meltdoy(matchtest), amj5cmSMmeanN,  '.g', ...
    meltdoy(matchtest), amj20cmSMmeanN, '.b', ...
    meltdoy(matchtest), amj50cmSMmeanN, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Day of wateryear');
ylabel('VWC (%)');
title(' Apr, May, Jun VWC vs snowmelt day');

subplot (2, 2, 4)
plot(meltdoy(matchtest), jas5cmSMmeanN,  '.g', ...
    meltdoy(matchtest), jas20cmSMmeanN, '.b', ...
    meltdoy(matchtest), jas50cmSMmeanN, '.k');
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
y = jas5cmSMmean;
y2 = jas5cmSMsd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(321); % 5cm  
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([140 310]); ylim([0 45]);
title('Summer VWC vs snowmelt day'); 
ylabel('5cm VWC');

y = jas20cmSMmean;
y2 = jas20cmSMsd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(323); % 20cm  
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([140 310]); ylim([0 45]);
ylabel('20cm VWC');

% 50cm
y = jas50cmSMmean;
y2 = jas50cmSMsd;
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
y = jas5cmSMmean;
y2 = jas5cmSMsd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(322);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([0 1650]); ylim([0 45]);
title('Summer VWC vs peak SWE');
ylabel('5cm VWC');

y = jas20cmSMmean;
y2 = jas20cmSMsd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(324);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([0 1650]); ylim([0 45]);
ylabel('20cm VWC');

y = jas50cmSMmean;
y2 = jas50cmSMsd;
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
y = jas5cmSMmeanN;
y2 = jas5cmSMsdN;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(321); % 5cm  
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([140 310]); ylim([0 1]);
title('Summer VWC vs snowmelt day'); 
ylabel('5cm VWC-norm');

y = jas20cmSMmeanN;
y2 = jas20cmSMsdN;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(323); % 20cm  
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([140 310]); ylim([0 1]);
ylabel('20cm VWC-norm');

% 50cm
y = jas50cmSMmeanN;
y2 = jas50cmSMsdN;
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
y = jas5cmSMmeanN;
y2 = jas5cmSMsdN;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(322);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([0 1650]); ylim([0 1]);
title('Summer VWC vs peak SWE'); 
ylabel('5cm VWC');

y = jas20cmSMmeanN;
y2 = jas20cmSMsdN;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(324);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([0 1650]); ylim([0 1]);
ylabel('20cm VWC');

y = jas50cmSMmeanN;
y2 = jas50cmSMsdN;
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
y = jas5cmSMmeanN;
subplot(321); % 5cm  
errorbar(x, y, jas5cmSMsdN, 'o', 'Color', [0.8 0.8 0.8],...
    'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', 'w');
hold on;
[~, rsq, xfit, yfit] = fitline(x, y, polyorder, [100, 400]);
plot(xfit, yfit, '-k', 'LineWidth', 2);
text(150, 0.75, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlim([140 310]); ylim([0 1]);
title('Summer VWC vs snowmelt day'); 
ylabel('5cm VWC-norm');

y = jas20cmSMmeanN;
subplot(323); % 20cm  
errorbar(x, y, jas20cmSMsdN, 'o', 'Color', [0.8 0.8 0.8],...
    'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', 'w');
hold on;
[~, rsq, xfit, yfit] = fitline(x, y, polyorder, [100, 400]);
plot(xfit, yfit, '-k', 'LineWidth', 2);
text(150, 0.75, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlim([140 310]); ylim([0 1]);
ylabel('20cm VWC-norm');

y = jas50cmSMmeanN;
subplot(325); % 50cm  
errorbar(x, y, jas50cmSMsdN, 'o', 'Color', [0.8 0.8 0.8],...
    'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', 'w');
hold on
[~, rsq, xfit, yfit] = fitline(x, y, polyorder, [100, 400]);
plot(xfit, yfit, '-k', 'LineWidth', 2);
text(150, 0.75, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlim([140 310]); ylim([0 1]);
ylabel('50cm VWC-norm');

% Right side - plot vs max SWE

x = maxswe(matchtest); %split into x and y
y = jas5cmSMmeanN;
subplot(322); % 5cm  
errorbar(x, y, jas5cmSMsdN, 'o', 'Color', [0.8 0.8 0.8],...
    'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', 'w');
hold on;
[~, rsq, xfit, yfit] = fitline(x, y, polyorder, [0, 1700]);
plot(xfit, yfit, '-k', 'LineWidth', 2);
text(150, 0.75, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlim([0 1650]); ylim([0 1]);
title('Summer VWC vs peak SWE'); 

y = jas20cmSMmeanN;
subplot(324); % 20cm  
errorbar(x, y, jas20cmSMsdN, 'ok', 'Color', [0.8 0.8 0.8],...
    'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', 'w');
hold on;
[~, rsq, xfit, yfit] = fitline(x, y, polyorder, [0, 1700]);
plot(xfit, yfit, '-k', 'LineWidth', 2);
text(150, 0.75, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlim([0 1650]); ylim([0 1]);

y = jas50cmSMmeanN;
subplot(326); % 50cm  
errorbar(x, y, jas50cmSMsdN, 'o', 'Color', [0.8 0.8 0.8],...
    'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', 'w');
hold on
[~, rsq, xfit, yfit] = fitline(x, y, polyorder, [0, 1700]);
plot(xfit, yfit, '-k', 'LineWidth', 2);
text(150, 0.75, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlim([0 1650]); ylim([0 1]);
xlabel('Peak SWE (mm)');

%--------------------------------------------------------------
% FIG 6 - Regress snowpack vs growing season vwc for 3 sites
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Regress vwc vs snowpack for 3 sites');
test = site_cl==828;
site1meltd = meltdoy(matchtest & test);
site1pks = maxswe(matchtest & test);
site1sm20 = jas20cmSMmeanN(site_sw==828);
test = site_cl==330;
site2meltd = meltdoy(matchtest & test);
site2pks = maxswe(matchtest & test);
site2sm20 = jas20cmSMmeanN(site_sw==330);
test = site_cl==582;
site3meltd = meltdoy(matchtest & test);
site3pks = maxswe(matchtest & test);
site3sm20 = jas20cmSMmeanN(site_sw==582);

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

%----------------------------------------------------
% FIG 7 - Plot August soil moisture vs snowmelt date
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Aug soil moisture vs snomelt date');

% Left side - plot vs Snowmelt day
% Set binning parameters
topEdge = 305; % define limits
botEdge = 145; % define limits
numBins = 15; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1))+5;

% 5cm
x = meltdoy(matchtest); %split into x and y
y = aug5cmSMmean;
y2 = aug5cmSMsd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(221);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([140 310]); %ylim([0 1]);
xlabel('Snowmelt dowy');
title('5cm VWC');

y = aug20cmSMmean;
y2 = aug20cmSMsd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(222);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([140 310]);%ylim([0 1]);
xlabel('Snowmelt dowy');
title('20cm VWC');

y = aug50cmSMmean;
y2 = aug50cmSMsd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins);
subplot(223);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([140 310]);%ylim([0 1]);
xlabel('Snowmelt dowy');
title('50cm VWC');

end