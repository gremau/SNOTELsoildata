function plot_snotel_summary
%plot_snotel_summary.m
%
% Reads the outputs from summarize_wateryear.m and makes a number of plots
% characterizing variability in soil moisture and soil temperature in the
% SNOTEL network.
%
% Feb 20, 2012 - Greg Maurer

clear;          % clear memory
close all;      % clear any figures
fignum=0;       % used to increment figure number for plots

% Ask user whether use normalized soil moisture data
normalize = input('Use normalized soil moisture data?  (y/n) : ', 's');

% Access to nan stuff, lines, etc
addpath('/home/greg/data/code_resources/m_common/nanstuff/');
addpath('/home/greg/data/code_resources/m_common/');
addpath('/home/greg/data/code_resources/m_common/hline_vline/'); 
%addpath('/home/greg/data/code_resources/m_common/hist2/');

% Set data path and file name, read in file
rawdatapath = '../rawdata/soilsensors_hourly/';
processeddatapath = '../processed_data/';

% Load lists of sites with data in the daily/hourly data directory
dailysites = sortrows(csvread('../rawdata/allsensors_daily/filelist.txt'));
soilsites = sortrows(csvread('../rawdata/soilsensors_hourly/filelist.txt'));
sites = unique(soilsites(:,1));

% LOAD the data (can switch between daily/hourly data here)
climData = csvread([processeddatapath 'wyear_climatesummary.txt']);
tsData = csvread([processeddatapath 'wyear_soiltempsummary_hourly.txt']);
% tsData = csvread([processeddatapath 'wyear_soiltempsummary_daily.txt']);
% Soil water content data
if strcmpi(normalize, 'y')
    vwcData = csvread([processeddatapath ...
        'wyear_soilwatersummary_hourly_smnorm.txt']);
%     vwcData = csvread([processeddatapath ...
%         'wyear_soilwatersummary_daily_smnorm.txt']);
else
    vwcData = csvread([processeddatapath...
        'wyear_soilwatersummary_hourly.txt']);
%     vwcData = csvread([processeddatapath...
%         'wyear_soilwatersummary_daily.txt']);
end

% climData includes more than just soil sites, 
% Get a subset that corresponds with available soildata
matchsoil = ismember(climData(:, 1:2), tsData(:, 1:2), 'rows');
soilClim = climData(matchsoil, :);


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

octTairMean = climData(:, 50);
octTairSd = climData(:, 51);
novTairMean = climData(:, 52);
novTairSd = climData(:, 53);
decTairMean = climData(:, 54);
decTairSd = climData(:, 55);
janTairMean = climData(:, 56);
janTairSd = climData(:, 57);
febTairMean = climData(:, 58);
febTairSd = climData(:, 59);
marTairMean = climData(:, 60);
marTairSd = climData(:, 61);
aprTairMean = climData(:, 62);
aprTairSd = climData(:, 63);
mayTairMean = climData(:, 64);
mayTairSd = climData(:, 65);
junTairMean = climData(:, 66);
junTairSd = climData(:, 67);
julTairMean = climData(:, 68);
julTairSd = climData(:, 69);
augTairMean = climData(:, 70);
augTairSd = climData(:, 71);
sepTairMean = climData(:, 72);
sepTairSd = climData(:, 73);
maat = climData(:, 74);
maat_sd = climData(:, 75);

preonsetTair = climData(:, 76);
preonsetTairSd = climData(:, 77);
premeltTair = climData(:, 78);
premeltTairSd = climData(:, 79);
postmeltTair = climData(:, 80);
postmeltTairSd = climData(:, 81);
elev = climData(:, 82);
lat = climData(:, 83);
lon = climData(:, 84);
ltMeanSWE = climData(:, 85);
ltMeanPrecip = climData(:, 86);

% Parse soilwatersummary
site_sw = vwcData(:, 1);
year_sw = vwcData(:, 2);
octVWC5mean = vwcData(:, 3);
octVWC5sd = vwcData(:, 4);
octVWC20mean = vwcData(:, 5);
octVWC20sd = vwcData(:, 6);
octVWC50mean = vwcData(:, 7);
octVWC50sd = vwcData(:, 8);
novVWC5mean = vwcData(:, 9);
novVWC5sd = vwcData(:, 10);
novVWC20mean = vwcData(:, 11);
novVWC20sd = vwcData(:, 12);
novVWC50mean = vwcData(:, 13);
novVWC50sd = vwcData(:, 14);
decVWC5mean = vwcData(:, 15);
decVWC5sd = vwcData(:, 16);
decVWC20mean = vwcData(:, 17);
decVWC20sd = vwcData(:, 18);
decVWC50mean = vwcData(:, 19);
decVWC50sd = vwcData(:, 20);
janVWC5mean = vwcData(:, 21);
janVWC5sd = vwcData(:, 22);
janVWC20mean = vwcData(:, 23);
janVWC20sd = vwcData(:, 24);
janVWC50mean = vwcData(:, 25);
janVWC50sd = vwcData(:, 26);
febVWC5mean = vwcData(:, 27);
febVWC5sd = vwcData(:, 28);
febVWC20mean = vwcData(:, 29);
febVWC20sd = vwcData(:, 30);
febVWC50mean = vwcData(:, 31);
febVWC50sd = vwcData(:, 32);
marVWC5mean = vwcData(:, 33);
marVWC5sd = vwcData(:, 34);
marVWC20mean = vwcData(:, 35);
marVWC20sd = vwcData(:, 36);
marVWC50mean = vwcData(:, 37);
marVWC50sd = vwcData(:, 38);
aprVWC5mean = vwcData(:, 39);
aprVWC5sd = vwcData(:, 40);
aprVWC20mean = vwcData(:, 41);
aprVWC20sd = vwcData(:, 42);
aprVWC50mean = vwcData(:, 43);
aprVWC50sd = vwcData(:, 44);

mayVWC5mean = vwcData(:, 45);
mayVWC5sd = vwcData(:, 46);
mayVWC20mean = vwcData(:, 47);
mayVWC20sd = vwcData(:, 48);
mayVWC50mean = vwcData(:, 49);
mayVWC50sd = vwcData(:, 50);

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

preonsetVWC5 = vwcData(:, 99);
preonsetVWC5sd = vwcData(:, 100);
preonsetVWC20 = vwcData(:, 101);
preonsetVWC20sd = vwcData(:, 102);
preonsetVWC50 = vwcData(:, 103);
preonsetVWC50sd = vwcData(:, 104);

% Parse soiltempsummary
site_st = tsData(:, 1);
year_st = tsData(:, 2);
octTs5mean = tsData(:, 3);
octTs5sd = tsData(:, 4);
octTs20mean = tsData(:, 5);
octTs20sd = tsData(:, 6);
octTs50mean = tsData(:, 7);
octTs50sd = tsData(:, 8);
decTs5mean = tsData(:, 15);
decTs5sd = tsData(:, 16);
decTs20mean = tsData(:, 17);
decTs20sd = tsData(:, 18);
decTs50mean = tsData(:, 19);
decTs50sd = tsData(:, 20);
% These repeat through sept (end of wy)
ondTs5mean = tsData(:, 75);
ondTs5sd = tsData(:, 76);
ondTs20mean = tsData(:, 77);
ondTs20sd = tsData(:, 78);
ondTs50mean = tsData(:, 79);
ondTs50sd = tsData(:, 80);
jfmTs5mean = tsData(:, 81);
jfmTs5sd = tsData(:, 82);
jfmTs20mean = tsData(:, 83);
jfmTs20sd = tsData(:, 84);
jfmTs50mean = tsData(:, 85);
jfmTs50sd = tsData(:, 86);
amjTs5mean = tsData(:, 87);
amjTs5sd = tsData(:, 88);
amjTs20mean = tsData(:, 89);
amjTs20sd = tsData(:, 90);
amjTs50mean = tsData(:, 91);
amjTs50sd = tsData(:, 92);
jasTs5mean = tsData(:, 93);
jasTs5sd = tsData(:, 94);
jasTs20mean = tsData(:, 95);
jasTs20sd = tsData(:, 96);
jasTs50mean = tsData(:, 97);
jasTs50sd = tsData(:, 98);

% Seasonal/yearly soil temp means
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

% Snowfree soil temp means
snowfreeTs5mean = tsData(:, 111);
snowfreeTs5sd = tsData(:, 112);
snowfreeTs20mean = tsData(:, 113);
snowfreeTs20sd = tsData(:, 114);
snowfreeTs50mean = tsData(:, 115);
snowfreeTs50sd = tsData(:, 116);

preonsetTs5 = tsData(:, 117);
preonsetTs5sd = tsData(:, 118);
preonsetTs20 = tsData(:, 119);
preonsetTs20sd = tsData(:, 120);
preonsetTs50 = tsData(:, 121);
preonsetTs50sd = tsData(:, 122);
premeltTs5 = tsData(:, 123);
premeltTs5sd = tsData(:, 124);
postmeltTs5 = tsData(:, 125);
postmeltTs5sd = tsData(:, 126);

%------------------------------------------------------------------
% FIG 1 - Plot data for entire network and soil subset

% Get the subset that is 
testsoil = ismember(site_cl, soilsites);

fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Wateryear data scatter - all sites/years');

subplot (4, 2, 1)
plot(elev, elev, 'ok', 'MarkerFaceColor', 'k');
hold on;
plot(elev(testsoil), elev(testsoil)+100, 'or');
ylabel('Elevation (m)');
legend('Intermountain west', 'Soil sites (+100m)');
title('Elevation');

subplot (4, 2, 2)
plot(elev, maat, 'ok', 'MarkerFaceColor', 'k');
hold on;
plot(elev(testsoil), maat(testsoil), 'or');
ylabel('MAT (deg C)');
legend('Intermountain west', 'Soil sites');
title('Mean wateryear air T');

subplot (4, 2, 3)
plot(elev, accumprecip, 'ok', 'MarkerFaceColor', 'k');
hold on;
testsoil = ismember(site_cl, soilsites);
plot(elev(testsoil), accumprecip(testsoil), 'or');
ylabel('Annual precip (mm)');
title('Wateryear precip');

subplot (4, 2, 4)
plot(elev, JASprecip, 'ok', 'MarkerFaceColor', 'k');
hold on;
plot(elev(testsoil), JASprecip(testsoil), 'or');
ylabel('Precip (mm)');
title('Summer Precip (Jul, Aug, Sep)');

subplot (4, 2, 5)
plot(elev, maxswe, 'ok', 'MarkerFaceColor', 'k');
hold on;
plot(elev(testsoil), maxswe(testsoil), 'or');
ylabel('SWE (mm)');
title('Peak SWE');

subplot (4, 2, 6)
plot(elev, totaldaysSC, 'ok', 'MarkerFaceColor', 'k');
hold on;
plot(elev(testsoil), totaldaysSC(testsoil), 'or');
ylabel('No. Days');
title('Total snowcovered days');

subplot (4, 2, 7)
plot(elev, onsetdoy, 'ok', 'MarkerFaceColor', 'k');
hold on;
plot(elev(testsoil), onsetdoy(testsoil), 'or');
xlabel('Elevation (m)'); ylabel('Day of year');
title('Snowpack onset day');

subplot (4, 2, 8)
plot(elev, meltdoy, 'ok', 'MarkerFaceColor', 'k');
hold on;
plot(elev(testsoil), meltdoy(testsoil), 'or');
xlabel('Elevation (m)'); ylabel('Day of year');
title('Day of snowmelt');

%------------------------------------------------------------------
% FIG 2 - Histograms of entire network and soil subset
testsoil = ismember(site_cl, soilsites);
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Wateryear data histograms - all sites/years');

subplot (4, 2, 1)
xedges = linspace(500, 4000, 60);
networkhist = histc(elev, xedges);
soilhist = histc(elev(testsoil), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'r');
vline(nanmean(elev), ':k');
vline(nanmean(elev(testsoil)), ':r');
xlim([700 3700]); ylim([0 400]);
title('Elevation');

subplot (4, 2, 2)
xedges = linspace(-5, 20, 60);
networkhist = histc(maat, xedges);
soilhist = histc(maat(testsoil), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'r');
vline(nanmean(maat), ':k');
vline(nanmean(maat(testsoil)), ':r');
xlim([-5 15]); ylim([0 500]);
title('Mean wateryear air T');

subplot (4, 2, 3)
xedges = linspace(0, 2500, 60);
networkhist = histc(accumprecip, xedges);
soilhist = histc(accumprecip(testsoil), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'r');
vline(nanmean(accumprecip), ':k');
vline(nanmean(accumprecip(testsoil)), ':r');
xlim([-10 2500]); ylim([0 500]);
title('Wateryear precip');

subplot (4, 2, 4)
xedges = linspace(0, 70, 60);
networkhist = histc(JASprecip, xedges);
soilhist = histc(JASprecip(testsoil), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'r');
vline(nanmean(JASprecip), ':k');
vline(nanmean(JASprecip(testsoil)), ':r');
xlim([-2 50]); ylim([0 700]);
title('Summer Precip (Jul, Aug, Sep)');

subplot (4, 2, 5)
xedges = linspace(100, 2000, 60);
networkhist = histc(maxswe, xedges);
soilhist = histc(maxswe(testsoil), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'r');
vline(nanmean(maxswe), ':k');
vline(nanmean(maxswe(testsoil)), ':r');
xlim([-2 2000]); ylim([0 400]);
title('Peak SWE');

subplot (4, 2, 6)
xedges = linspace(0, 365, 60);
networkhist = histc(totaldaysSC, xedges);
soilhist = histc(totaldaysSC(testsoil), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'r');
vline(nanmean(totaldaysSC), ':k');
vline(nanmean(totaldaysSC(testsoil)), ':r');
xlim([-2 300]); ylim([0 500]);
title('Total snowcovered days');

subplot (4, 2, 7)
xedges = linspace(0, 130, 60);
networkhist = histc(onsetdoy, xedges);
soilhist = histc(onsetdoy(testsoil), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'r');
vline(nanmean(onsetdoy), ':k');
vline(nanmean(onsetdoy(testsoil)), ':r');
xlim([-2 125]); ylim([0 700]);
title('Snowpack onset day');

subplot (4, 2, 8)
xedges = linspace(0, 365, 60);
networkhist = histc(meltdoy, xedges);
soilhist = histc(meltdoy(testsoil), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'r');
vline(nanmean(meltdoy), ':k');
vline(nanmean(meltdoy(testsoil)), ':r');
xlim([100 350]); ylim([0 700]);
title('Day of snowmelt');


% %------------------------------------------------------------------
% % FIG 3 - 2d histogram of some data above (just as an example)
% fignum = fignum+1;
% h = figure(fignum);
% set(h, 'Name', 'Wateryear snowpack metrics 1 - all sites');
% 
% colormap(jet);
% 
% subplot (2, 1, 1)
% % 2D histogram
% test = ~isnan(elev) & ~isnan(maat);
% x = elev(test);
% y = maat(test);
% xedges = linspace(500, 4000, 60);
% yedges = linspace(-2, 14, 60);
% histmat = hist2(x, y, xedges, yedges);  % hist2 is from the matlab user forum
% pcolor(xedges, yedges, histmat');
% title('Mean wateryear Tair across network');
% 
% subplot (2, 1, 2)
% plot(elev, maat, 'om', 'MarkerFaceColor', 'm');
% hold on;
% testsoil = ismember(site_cl, soilsites);
% plot(elev(testsoil), maat(testsoil), 'ok');
% xlabel('Elevation (m)');
% ylabel('Mean wateryear Tair (^oC)');
% legend('Intermountain west', 'Soil sites');

%--------------------------------------------------------------
% FIG 3 - Plot MAST & MAT vs max swe, snowcover duration, elevation
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Compare MAST & MAT (vs maxSWE, snowduration, elev');

% First plot MAT vs SWE, and its fit line/r-squared value
subplot1 = subplot (2, 2, 1);
plot(maxswe(matchsoil), maat(matchsoil), 'ok', ...
    'MarkerFaceColor', 'red');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(maxswe(matchsoil), maat(matchsoil), 2, xrange);
plot(xfit, yfit, '--k');
text(1200, 0, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
% Then plot MAST and its fit line/r-squared value
plot(maxswe(matchsoil), mast5cm, '.b', ...
    maxswe(matchsoil), mast20cm, '+b', ...
    maxswe(matchsoil), mast50cm, '*b');
[~, rsq, xfit, yfit] = fitline(maxswe(matchsoil), mast20cm, 2, xrange);
plot(xfit, yfit, ':k');
text(1200, 6, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
% Label stuff
legend('Mean Air T', 'linear fit', 'Mean 5cm Ts','20cm','50cm',...
    'linear fit (20cm)');
xlabel('Peak SWE (mm)'); ylabel('^oC');
title('Mean wateryear Tair & SoilT vs peak SWE');

% Then plot MAT vs elevation, and its fit line/r-squared value
subplot (2, 2, 2)
plot(elev(matchsoil), maat(matchsoil), 'ok', ...
    'MarkerFaceColor', 'red');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(elev(matchsoil), maat(matchsoil), 1, xrange);
plot(xfit, yfit, '--k');
text(500, 12, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
% MAST and its fit line/r-squared value
plot(elev(matchsoil), mast5cm, '.b', ...
    elev(matchsoil), mast20cm, '+b',...
    elev(matchsoil), mast50cm, '*b');
[~, rsq, xfit, yfit] = fitline(elev(matchsoil), mast20cm, 1, xrange);
plot(xfit, yfit, ':k');
text(500, 7, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
% Label stuff
xlabel('Elevation (m)'); ylabel('^oC');
title('Mean wateryear Tair & SoilT vs Elevation');

% Then MAT vs snowpack duration, and its fit line/r-squared value
subplot (2, 2, 3)
plot(snowduration(matchsoil), maat(matchsoil), 'ok', ...
    'MarkerFaceColor', 'red');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(snowduration(matchsoil), maat(matchsoil), 1, xrange);
plot(xfit, yfit, '--k');
text(30, 7, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
% MAST and its fit line/r-squared value
plot(snowduration(matchsoil), mast5cm, '.b', ...
    snowduration(matchsoil), mast20cm, '+b',...
    snowduration(matchsoil), mast50cm, '*b');
[~, rsq, xfit, yfit] = fitline(snowduration(matchsoil), mast20cm, 1, xrange);
plot(xfit, yfit, ':k');
text(30, 12, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
% Label stuff
xlabel('Snowpack duration (days)'); ylabel('^oC');
title('... vs Snowpack duration');

% And MAT vs total snowcovered days, and its fit line/r-squared value
subplot (2, 2, 4)
plot(totaldaysSC(matchsoil), maat(matchsoil), 'ok', ...
    'MarkerFaceColor', 'red');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(totaldaysSC(matchsoil), maat(matchsoil), 1, xrange);
plot(xfit, yfit, '--k');
text(30, 7, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
% MAST and its fit line/r-squared value
plot(totaldaysSC(matchsoil), mast5cm, '.b', ...
    totaldaysSC(matchsoil), mast20cm, '+b',...
    totaldaysSC(matchsoil), mast50cm, '*b');
[~, rsq, xfit, yfit] = fitline(totaldaysSC(matchsoil), mast20cm, 1, xrange);
plot(xfit, yfit, ':k');
text(30, 12, ['r^2 = ' num2str(rsq, 2)]);
% Label stuff
xlabel('No. snowcovered days'); ylabel('^oC');
title('... vs Tot. snowcovered days');

%--------------------------------------------------------------
% FIG 4 - Plot MAST vs MAT/Snowcovered days by depth and elevation
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'MAST vs MAT/Snowcovered days by depth and elevation');

% Left side - Plot MAST vs MAT/snowcovered days by depth
subplot (2, 2, 1)
plot(maat(matchsoil), mast5cm, '.g', ...
    maat(matchsoil), mast20cm, '.b',...
    maat(matchsoil), mast50cm, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Mean Tair (^oC)');
ylabel('MAST(^oC)');
title('Mean wateryear Ts vs Mean wateryear Tair');

subplot (2, 2, 3)
plot(totaldaysSC(matchsoil), mast5cm, '.g', ...
    totaldaysSC(matchsoil), mast20cm, '.b',...
    totaldaysSC(matchsoil), mast50cm, '.k');
xlabel('No. days');
ylabel('MAST (^oC)');
title('... vs Snowcovered days');

% Right side - Plot MAST vs MAT/Snowcovered days in elevation bins
% First get the elevation categories
testhi = elev(matchsoil) > 3000;
testmid = elev(matchsoil) > 2500 & elev(matchsoil) < 3000;
testlo = elev(matchsoil) > 2000 & elev(matchsoil) < 2500;
% Match the climatedata variables with soil data
matchSCdays = totaldaysSC(matchsoil);
matchMAT = maat(matchsoil);

subplot (2, 2, 2)
plot(matchMAT(testhi), mast20cm(testhi), '.b', ...
    matchMAT(testmid), mast20cm(testmid), '.g',...
    matchMAT(testlo), mast20cm(testlo), '.r');
legend('3000+', '2500-3000', '2000-2500cm');
xlabel('Mean Tair (^oC)');
ylabel('20cm MAST (^oC)');
title('Mean wateryear Ts vs. Mean wateryear Tair');

subplot (2, 2, 4)
plot(matchSCdays(testhi), mast20cm(testhi), '.b', ...
    matchSCdays(testmid), mast20cm(testmid), '.g',...
    matchSCdays(testlo), mast20cm(testlo), '.r');
xlabel('No. days');
ylabel('20cm MAST(^oC)');
title('... vs Snowcovered days');

%--------------------------------------------------------------
% FIG 5 - Plot OFFSETS between MAST and MAT
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'MAST-MAT (Temp offset) for all wateryears');
% Calculate the offsets for each depth
offset5cm = mast5cm-maat(matchsoil);
offset20cm = mast20cm-maat(matchsoil);
offset50cm = mast50cm-maat(matchsoil);

subplot (2, 2, 1)
plot(elev(matchsoil), offset5cm, '.r', ...
    elev(matchsoil), offset20cm, '.g',...
    elev(matchsoil), offset50cm, '.b');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(elev(matchsoil), offset20cm, 1, xrange);
plot(xfit, yfit, '--k');
text(1000, 7, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
legend('5cm', '20cm', '50cm', 'linear fit (20cm)');
xlabel('Elevation(m)'); ylabel('^oC');
title('MAST - MAT vs Elevation');

subplot (2, 2, 2)
plot(totaldaysSC(matchsoil), offset5cm, '.r', ...
    totaldaysSC(matchsoil), offset20cm, '.g',...
    totaldaysSC(matchsoil), offset50cm, '.b');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(totaldaysSC(matchsoil), offset20cm, 1, xrange);
plot(xfit, yfit, '--k');
text(30, 7, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('No. days'); ylabel('^oC');
title('MAST - MAT vs Tot. Snowcovered Days');

subplot (2, 2, 3)
plot(maxswe(matchsoil), offset5cm, '.r', ...
    maxswe(matchsoil), offset20cm, '.g',...
    maxswe(matchsoil), offset50cm, '.b');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(maxswe(matchsoil), offset20cm, 1, xrange);
plot(xfit, yfit, '--k');
text(200, 7, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('mm'); ylabel('^oC');
title('MAST - MAT vs Peak SWE');

%--------------------------------------------------------------
% FIG 6 - Plot MAST vs snowpack duration for 3 sites
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Regress mast vs snowpack duration for 3 sites');
test = site_cl==828;
site1dur = snowduration(matchsoil & test);
site1maat = maat(matchsoil & test);
site1mast = mast20cm(site_st==828);
test = site_cl==393;
site2dur = snowduration(matchsoil & test);
site2maat = maat(matchsoil & test);
site2mast = mast20cm(site_st==393);
test = site_cl==582;
site3dur = snowduration(matchsoil & test);
site3mast = mast20cm(site_st==582);
site3maat = maat(matchsoil & test);

subplot (3, 2, 1)
plot(site1dur, site1mast, '.b');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site1dur, site1mast, 1, xrange);
plot(xfit, yfit,':k');
text(160, 4, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('Snowpack duration (days)'); ylabel('MAST');
title('Trial Lake');

subplot (3, 2, 3)
plot(site2dur, site2mast, '.b');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site2dur, site2mast, 1, xrange);
plot(xfit, yfit,':k');
text(190, 5, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('Snowpack duration (days)'); ylabel('MAST');
title('Chalk Creek 1');

subplot (3, 2, 5)
plot(site3dur, site3mast, '.b');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site3dur, site3mast, 1, xrange);
plot(xfit, yfit,':k');
text(165, 8.5, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('Snowpack duration (days)'); ylabel('MAST');
title('Little Bear');

subplot (3, 2, 2)
plot(site1maat, site1mast, '.b')
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site1maat, site1mast, 1, xrange);
plot(xfit, yfit,':k');
text(1, 2.7, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('MAT'); ylabel('MAST');
title('Trial Lake');

subplot (3, 2, 4)
plot(site2maat, site2mast, '.b')
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site2maat, site2mast, 1, xrange);
plot(xfit, yfit,':k');
text(3.7, 4.5, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('MAT'); ylabel('MAST');
title('Chalk Creek 1');

subplot (3, 2, 6)
plot(site3maat, site3mast, '.b')
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site3maat, site3mast, 1, xrange);
plot(xfit, yfit,':k');
text(7, 7.5, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('MAT'); ylabel('MAST');
title('Little Bear');

%----------------------------------------------------
% FIG 7 - Plot snowcovered temp vs onset temps
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Winter SoilT vs pre-snowpack temps');

subplot (2, 2, 1)
plot(preonsetTair(matchsoil), snowcovTs5mean, '.g', ...
    preonsetTair(matchsoil), snowcovTs20mean, '.b', ...
    preonsetTair(matchsoil), snowcovTs50mean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Wateryear pre-snowpack Tair (^oC)');
ylabel('Mean wateryear belowsnow SoilT (^oC)');
title('Mean wateryear belowsnow SoilT vs. pre-snowpack Air T');

subplot (2, 2, 2)
plot(preonsetTs5, snowcovTs5mean, '.g', ...
    preonsetTs20, snowcovTs20mean, '.b', ...
    preonsetTs50, snowcovTs50mean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Wateryear pre-snowpack SoilT (^oC)');
ylabel('Mean wateryear belowsnow SoilT (^oC)');
title('Mean wateryear belowsnow SoilT vs. pre-snowpack SoilT');

subplot (2, 2, 3)
plot(preonsetTs5, ondTs5mean, '.g', ...
    preonsetTs20, ondTs20mean, '.b', ...
    preonsetTs50, ondTs50mean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Wateryear pre-snowpack SoilT (^oC)');
ylabel('Mean OND temp (^oC)');
title('Mean Oct, Nov, Dec SoilT vs. pre-snowpack SoilT');

subplot (2, 2, 4)
plot(preonsetTs5, jfmTs5mean, '.g', ...
    preonsetTs20, jfmTs20mean, '.b', ...
    preonsetTs50, jfmTs50mean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Wateryear pre-snowpack SoilT (^oC)');
ylabel('Mean JFM temp (^oC)');
title('Mean Jan, Feb, Mar SoilT vs. pre-snowpack SoilT');

%------------------------------------------------------------------
% FIG 8 - Plot soil moisture vs max swe/snowmelt day
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Wateryear seasonal soil moisture');

subplot (2, 2, 1)
plot(maxswe(matchsoil), amjVWC5mean,  '.g', ...
    maxswe(matchsoil), amjVWC20mean, '.b', ...
    maxswe(matchsoil), amjVWC50mean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Peak SWE (mm)');
ylabel('VWC (%)');
title(' Apr, May, Jun VWC');

subplot (2, 2, 2)
plot(maxswe(matchsoil), jasVWC5mean,  '.g', ...
    maxswe(matchsoil), jasVWC20mean, '.b', ...
    maxswe(matchsoil), jasVWC50mean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Peak SWE (mm)');
ylabel('VWC (%)');
title(' July, Aug, Sept VWC');

subplot (2, 2, 3)
plot(meltdoy(matchsoil), amjVWC5mean,  '.g', ...
    meltdoy(matchsoil), amjVWC20mean, '.b', ...
    meltdoy(matchsoil), amjVWC50mean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Day of wateryear');
ylabel('VWC (%)');
title(' Apr, May, Jun VWC vs snowmelt day');

subplot (2, 2, 4)
plot(meltdoy(matchsoil), jasVWC5mean,  '.g', ...
    meltdoy(matchsoil), jasVWC20mean, '.b', ...
    meltdoy(matchsoil), jasVWC50mean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Day of wateryear (mm)');
ylabel('VWC (%)');
title(' July, Aug, Sept VWC vs snowmeltday');

%----------------------------------------------------
% FIG 9 - Plot winter soil moisture vs onset soil moisture
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Winter VWC vs pre-snowpack VWC');

subplot (2, 2, 1)
plot(preonsetVWC5, decVWC5mean, '.g', ...
    preonsetVWC20, decVWC20mean, '.b', ...
    preonsetVWC50, decVWC50mean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Wateryear pre-snowpack soilVWC (%)');
ylabel('Mean VWC (%)');
title('Oct, Nov, Dec VWC vs pre-snowpack VWC');

subplot (2, 2, 2)
plot(preonsetVWC5, febVWC5mean, '.g', ...
    preonsetVWC20, febVWC20mean, '.b', ...
    preonsetVWC50, febVWC50mean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Wateryear pre-snowpack soilVWC (%)');
ylabel('Mean VWC (%)');
title('February VWC vs pre-snowpack VWC');

%----------------------------------------------------
% FIG 10 - Plot winter soil moisture vs onset soil moisture
% Same as above, but binned, using only 20cm data
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Winter VWC vs pre-snowpack VWC');


% Set binning parameters
if strcmpi(normalize, 'n')
    topEdge = 45; % define limits
    botEdge = 0; % define limits
    numBins = 10; % define number of bins
    xaxlim = [0 45];
    yaxlim = [0 45];
elseif strcmpi(normalize, 'y')
    topEdge = 1; % define limits
    botEdge = 0; % define limits
    numBins = 10; % define number of bins
    xaxlim = [-0.1 1];
    yaxlim = [0 1];
end
    
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1));

% And months to plot
months = ['Nov'; 'Jan'; 'Mar'];

for i = 1:3;
    subplot (3,1, i)
    x = preonsetVWC5; %split into x and y
    eval(['y = ' lower(months(i,:)) 'VWC5mean;']);
    eval(['y2 = ' lower(months(i,:)) 'VWC5sd;']);
    [binMeanDec, binSdDec] = binseries(x, y, y2, topEdge, botEdge, numBins);
    errorbar(xax(1:numBins), binMeanDec, binSdDec, 'o',...
        'Color', [0.7 0.7 0.7]);
    hold on;
    x = preonsetVWC20; %split into x and y
    eval(['y = ' lower(months(i,:)) 'VWC20mean;']);
    eval(['y2 = ' lower(months(i,:)) 'VWC20sd;']);
    [binMeanDec, binSdDec] = binseries(x, y, y2, topEdge, botEdge, numBins);
    errorbar(xax(1:numBins), binMeanDec, binSdDec, 'o',...
        'Color', [0.4 0.4 0.4]);
    x = preonsetVWC50; %split into x and y
    eval(['y = ' lower(months(i,:)) 'VWC50mean;']);
    eval(['y2 = ' lower(months(i,:)) 'VWC50sd;']);
    [binMeanDec, binSdDec] = binseries(x, y, y2, topEdge, botEdge, numBins);
    errorbar(xax(1:numBins), binMeanDec, binSdDec, 'ok');
    xlabel('Wateryear pre-snowpack soilVWC (%)');
    ylabel('Mean VWC (%)');
    xlim(xaxlim); ylim(yaxlim);
    title([months(i,:) ' VWC vs pre-snowpack VWC']);
    if i==1
        legend('5cm', '20cm', '50cm', 'Location', 'NorthWest');
    end
    
end
junk = 99;
end



