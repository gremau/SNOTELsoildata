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

% LOAD the data (can switch between daily/hourly data here
climData = csvread([processeddatapath 'wyear_climatesummary.txt']);
filelistSoil = dlmread([rawdatapath 'filelist.txt'], ',');
soilsites = unique(filelistSoil(:, 1));
% Soil temp data
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
nov5cmSMmean = vwcData(:, 9);
nov5cmSMsd = vwcData(:, 10);
nov20cmSMmean = vwcData(:, 11);
nov20cmSMsd = vwcData(:, 12);
nov50cmSMmean = vwcData(:, 13);
nov50cmSMsd = vwcData(:, 14);
dec5cmSMmean = vwcData(:, 15);
dec5cmSMsd = vwcData(:, 16);
dec20cmSMmean = vwcData(:, 17);
dec20cmSMsd = vwcData(:, 18);
dec50cmSMmean = vwcData(:, 19);
dec50cmSMsd = vwcData(:, 20);
jan5cmSMmean = vwcData(:, 21);
jan5cmSMsd = vwcData(:, 22);
jan20cmSMmean = vwcData(:, 23);
jan20cmSMsd = vwcData(:, 24);
jan50cmSMmean = vwcData(:, 25);
jan50cmSMsd = vwcData(:, 26);
feb5cmSMmean = vwcData(:, 27);
feb5cmSMsd = vwcData(:, 28);
feb20cmSMmean = vwcData(:, 29);
feb20cmSMsd = vwcData(:, 30);
feb50cmSMmean = vwcData(:, 31);
feb50cmSMsd = vwcData(:, 32);
mar5cmSMmean = vwcData(:, 33);
mar5cmSMsd = vwcData(:, 34);
mar20cmSMmean = vwcData(:, 35);
mar20cmSMsd = vwcData(:, 36);
mar50cmSMmean = vwcData(:, 37);
mar50cmSMsd = vwcData(:, 38);
apr5cmSMmean = vwcData(:, 39);
apr5cmSMsd = vwcData(:, 40);
apr20cmSMmean = vwcData(:, 41);
apr20cmSMsd = vwcData(:, 42);
apr50cmSMmean = vwcData(:, 43);
apr50cmSMsd = vwcData(:, 44);

may5cmSMmean = vwcData(:, 45);
may5cmSMsd = vwcData(:, 46);
may20cmSMmean = vwcData(:, 47);
may20cmSMsd = vwcData(:, 48);
may50cmSMmean = vwcData(:, 49);
may50cmSMsd = vwcData(:, 50);

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

preonset5cmSM = vwcData(:, 99);
preonset5cmSMsd = vwcData(:, 100);
preonset20cmSM = vwcData(:, 101);
preonset20cmSMsd = vwcData(:, 102);
preonset50cmSM = vwcData(:, 103);
preonset50cmSMsd = vwcData(:, 104);

% Parse soiltempsummary
site_st = tsData(:, 1);
year_st = tsData(:, 2);
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
ond5cmSTmean = tsData(:, 75);
ond5cmSTsd = tsData(:, 76);
ond20cmSTmean = tsData(:, 77);
ond20cmSTsd = tsData(:, 78);
ond50cmSTmean = tsData(:, 79);
ond50cmSTsd = tsData(:, 80);
jfm5cmSTmean = tsData(:, 81);
jfm5cmSTsd = tsData(:, 82);
jfm20cmSTmean = tsData(:, 83);
jfm20cmSTsd = tsData(:, 84);
jfm50cmSTmean = tsData(:, 85);
jfm50cmSTsd = tsData(:, 86);
amj5cmSTmean = tsData(:, 87);
amj5cmSTsd = tsData(:, 88);
amj20cmSTmean = tsData(:, 89);
amj20cmSTsd = tsData(:, 90);
amj50cmSTmean = tsData(:, 91);
amj50cmSTsd = tsData(:, 92);
jas5cmSTmean = tsData(:, 93);
jas5cmSTsd = tsData(:, 94);
jas20cmSTmean = tsData(:, 95);
jas20cmSTsd = tsData(:, 96);
jas50cmSTmean = tsData(:, 97);
jas50cmSTsd = tsData(:, 98);

% Seasonal/yearly soil temp means
mast5cm = tsData(:, 99);
sdast5cm = tsData(:, 100);
mast20cm = tsData(:, 101);
sdast20cm = tsData(:, 102);
mast50cm = tsData(:, 103);
sdast50cm = tsData(:, 104);

% Snowcovered soil temp means
snowcovMeanST5cm = tsData(:, 105);
snowcovStdST5cm = tsData(:, 106);
snowcovMeanST20cm = tsData(:, 107);
snowcovStdST20cm = tsData(:, 108);
snowcovMeanST50cm = tsData(:, 109);
snowcovStdST50cm = tsData(:, 110);

% Snowfree soil temp means
snowfreeMeanST5cm = tsData(:, 111);
snowfreeStdST5cm = tsData(:, 112);
snowfreeMeanST20cm = tsData(:, 113);
snowfreeStdST20cm = tsData(:, 114);
snowfreeMeanST50cm = tsData(:, 115);
snowfreeStdST50cm = tsData(:, 116);

preonset5cmST = tsData(:, 117);
preonset5cmSTsd = tsData(:, 118);
preonset20cmST = tsData(:, 119);
preonset20cmSTsd = tsData(:, 120);
preonset50cmST = tsData(:, 121);
preonset50cmSTsd = tsData(:, 122);
premelt5cmST = tsData(:, 123);
premelt5cmSTsd = tsData(:, 124);
postmelt5cmST = tsData(:, 125);
postmelt5cmSTsd = tsData(:, 126);

% Get a subset of climData that corresponds with available soildata
matchtest = ismember(climData(:, 1:2), filelistSoil(:, 1:2), 'rows');
matchsets = climData(matchtest, :);

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
% title('Mean wateryear airT across network');
% 
% subplot (2, 1, 2)
% plot(elev, maat, 'om', 'MarkerFaceColor', 'm');
% hold on;
% testsoil = ismember(site_cl, soilsites);
% plot(elev(testsoil), maat(testsoil), 'ok');
% xlabel('Elevation (m)');
% ylabel('Mean wateryear airT (^oC)');
% legend('Intermountain west', 'Soil sites');

%--------------------------------------------------------------
% FIG 3 - Plot MAST & MAT vs max swe, snowcover duration, elevation
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Compare MAST & MAT (vs maxSWE, snowduration, elev');

% First plot MAT vs SWE, and its fit line/r-squared value
subplot1 = subplot (2, 2, 1);
plot(maxswe(matchtest), maat(matchtest), 'ok', ...
    'MarkerFaceColor', 'red');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(maxswe(matchtest), maat(matchtest), 2, xrange);
plot(xfit, yfit, '--k');
text(1200, 0, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
% Then plot MAST and its fit line/r-squared value
plot(maxswe(matchtest), mast5cm, '.b', ...
    maxswe(matchtest), mast20cm, '+b', ...
    maxswe(matchtest), mast50cm, '*b');
[~, rsq, xfit, yfit] = fitline(maxswe(matchtest), mast20cm, 2, xrange);
plot(xfit, yfit, ':k');
text(1200, 6, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
% Label stuff
legend('Mean Air T', 'linear fit', 'Mean 5cm Ts','20cm','50cm',...
    'linear fit (20cm)');
xlabel('Peak SWE (mm)'); ylabel('^oC');
title('Mean wateryear AirT & SoilT vs peak SWE');

% Then plot MAT vs elevation, and its fit line/r-squared value
subplot (2, 2, 2)
plot(elev(matchtest), maat(matchtest), 'ok', ...
    'MarkerFaceColor', 'red');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(elev(matchtest), maat(matchtest), 1, xrange);
plot(xfit, yfit, '--k');
text(500, 12, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
% MAST and its fit line/r-squared value
plot(elev(matchtest), mast5cm, '.b', ...
    elev(matchtest), mast20cm, '+b',...
    elev(matchtest), mast50cm, '*b');
[~, rsq, xfit, yfit] = fitline(elev(matchtest), mast20cm, 1, xrange);
plot(xfit, yfit, ':k');
text(500, 7, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
% Label stuff
xlabel('Elevation (m)'); ylabel('^oC');
title('Mean wateryear AirT & SoilT vs Elevation');

% Then MAT vs snowpack duration, and its fit line/r-squared value
subplot (2, 2, 3)
plot(snowduration(matchtest), maat(matchtest), 'ok', ...
    'MarkerFaceColor', 'red');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(snowduration(matchtest), maat(matchtest), 1, xrange);
plot(xfit, yfit, '--k');
text(30, 7, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
% MAST and its fit line/r-squared value
plot(snowduration(matchtest), mast5cm, '.b', ...
    snowduration(matchtest), mast20cm, '+b',...
    snowduration(matchtest), mast50cm, '*b');
[~, rsq, xfit, yfit] = fitline(snowduration(matchtest), mast20cm, 1, xrange);
plot(xfit, yfit, ':k');
text(30, 12, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
% Label stuff
xlabel('Snowpack duration (days)'); ylabel('^oC');
title('... vs Snowpack duration');

% And MAT vs total snowcovered days, and its fit line/r-squared value
subplot (2, 2, 4)
plot(totaldaysSC(matchtest), maat(matchtest), 'ok', ...
    'MarkerFaceColor', 'red');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(totaldaysSC(matchtest), maat(matchtest), 1, xrange);
plot(xfit, yfit, '--k');
text(30, 7, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
% MAST and its fit line/r-squared value
plot(totaldaysSC(matchtest), mast5cm, '.b', ...
    totaldaysSC(matchtest), mast20cm, '+b',...
    totaldaysSC(matchtest), mast50cm, '*b');
[~, rsq, xfit, yfit] = fitline(totaldaysSC(matchtest), mast20cm, 1, xrange);
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
plot(maat(matchtest), mast5cm, '.g', ...
    maat(matchtest), mast20cm, '.b',...
    maat(matchtest), mast50cm, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Mean AirT (^oC)');
ylabel('MAST(^oC)');
title('Mean wateryear Ts vs Mean wateryear Tair');

subplot (2, 2, 3)
plot(totaldaysSC(matchtest), mast5cm, '.g', ...
    totaldaysSC(matchtest), mast20cm, '.b',...
    totaldaysSC(matchtest), mast50cm, '.k');
xlabel('No. days');
ylabel('MAST (^oC)');
title('... vs Snowcovered days');

% Right side - Plot MAST vs MAT/Snowcovered days in elevation bins
% First get the elevation categories
testhi = elev(matchtest) > 3000;
testmid = elev(matchtest) > 2500 & elev(matchtest) < 3000;
testlo = elev(matchtest) > 2000 & elev(matchtest) < 2500;
% Match the climatedata variables with soil data
matchSCdays = totaldaysSC(matchtest);
matchMAT = maat(matchtest);

subplot (2, 2, 2)
plot(matchMAT(testhi), mast20cm(testhi), '.b', ...
    matchMAT(testmid), mast20cm(testmid), '.g',...
    matchMAT(testlo), mast20cm(testlo), '.r');
legend('3000+', '2500-3000', '2000-2500cm');
xlabel('Mean AirT (^oC)');
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
offset5cm = mast5cm-maat(matchtest);
offset20cm = mast20cm-maat(matchtest);
offset50cm = mast50cm-maat(matchtest);

subplot (2, 2, 1)
plot(elev(matchtest), offset5cm, '.r', ...
    elev(matchtest), offset20cm, '.g',...
    elev(matchtest), offset50cm, '.b');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(elev(matchtest), offset20cm, 1, xrange);
plot(xfit, yfit, '--k');
text(1000, 7, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
legend('5cm', '20cm', '50cm', 'linear fit (20cm)');
xlabel('Elevation(m)'); ylabel('^oC');
title('MAST - MAT vs Elevation');

subplot (2, 2, 2)
plot(totaldaysSC(matchtest), offset5cm, '.r', ...
    totaldaysSC(matchtest), offset20cm, '.g',...
    totaldaysSC(matchtest), offset50cm, '.b');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(totaldaysSC(matchtest), offset20cm, 1, xrange);
plot(xfit, yfit, '--k');
text(30, 7, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('No. days'); ylabel('^oC');
title('MAST - MAT vs Tot. Snowcovered Days');

subplot (2, 2, 3)
plot(maxswe(matchtest), offset5cm, '.r', ...
    maxswe(matchtest), offset20cm, '.g',...
    maxswe(matchtest), offset50cm, '.b');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(maxswe(matchtest), offset20cm, 1, xrange);
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
site1dur = snowduration(matchtest & test);
site1maat = maat(matchtest & test);
site1mast = mast20cm(site_st==828);
test = site_cl==393;
site2dur = snowduration(matchtest & test);
site2maat = maat(matchtest & test);
site2mast = mast20cm(site_st==393);
test = site_cl==582;
site3dur = snowduration(matchtest & test);
site3mast = mast20cm(site_st==582);
site3maat = maat(matchtest & test);

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
plot(preonsetAirT(matchtest), snowcovMeanST5cm, '.g', ...
    preonsetAirT(matchtest), snowcovMeanST20cm, '.b', ...
    preonsetAirT(matchtest), snowcovMeanST50cm, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Wateryear pre-snowpack AirT (^oC)');
ylabel('Mean wateryear belowsnow SoilT (^oC)');
title('Mean wateryear belowsnow SoilT vs. pre-snowpack Air T');

subplot (2, 2, 2)
plot(preonset5cmST, snowcovMeanST5cm, '.g', ...
    preonset20cmST, snowcovMeanST20cm, '.b', ...
    preonset50cmST, snowcovMeanST50cm, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Wateryear pre-snowpack SoilT (^oC)');
ylabel('Mean wateryear belowsnow SoilT (^oC)');
title('Mean wateryear belowsnow SoilT vs. pre-snowpack SoilT');

subplot (2, 2, 3)
plot(preonset5cmST, ond5cmSTmean, '.g', ...
    preonset20cmST, ond20cmSTmean, '.b', ...
    preonset50cmST, ond50cmSTmean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Wateryear pre-snowpack SoilT (^oC)');
ylabel('Mean OND temp (^oC)');
title('Mean Oct, Nov, Dec SoilT vs. pre-snowpack SoilT');

subplot (2, 2, 4)
plot(preonset5cmST, jfm5cmSTmean, '.g', ...
    preonset20cmST, jfm20cmSTmean, '.b', ...
    preonset50cmST, jfm50cmSTmean, '.k');
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
% FIG 9 - Plot winter soil moisture vs onset soil moisture
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Winter VWC vs pre-snowpack VWC');

subplot (2, 2, 1)
plot(preonset5cmSM, dec5cmSMmean, '.g', ...
    preonset20cmSM, dec20cmSMmean, '.b', ...
    preonset50cmSM, dec50cmSMmean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Wateryear pre-snowpack soilVWC (%)');
ylabel('Mean VWC (%)');
title('Oct, Nov, Dec VWC vs pre-snowpack VWC');

subplot (2, 2, 2)
plot(preonset5cmSM, feb5cmSMmean, '.g', ...
    preonset20cmSM, feb20cmSMmean, '.b', ...
    preonset50cmSM, feb50cmSMmean, '.k');
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
    x = preonset5cmSM; %split into x and y
    eval(['y = ' lower(months(i,:)) '5cmSMmean;']);
    eval(['y2 = ' lower(months(i,:)) '5cmSMsd;']);
    [binMeanDec, binSdDec] = binseries(x, y, y2, topEdge, botEdge, numBins);
    errorbar(xax(1:numBins), binMeanDec, binSdDec, 'o',...
        'Color', [0.7 0.7 0.7]);
    hold on;
    x = preonset20cmSM; %split into x and y
    eval(['y = ' lower(months(i,:)) '20cmSMmean;']);
    eval(['y2 = ' lower(months(i,:)) '20cmSMsd;']);
    [binMeanDec, binSdDec] = binseries(x, y, y2, topEdge, botEdge, numBins);
    errorbar(xax(1:numBins), binMeanDec, binSdDec, 'o',...
        'Color', [0.4 0.4 0.4]);
    x = preonset50cmSM; %split into x and y
    eval(['y = ' lower(months(i,:)) '50cmSMmean;']);
    eval(['y2 = ' lower(months(i,:)) '50cmSMsd;']);
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



