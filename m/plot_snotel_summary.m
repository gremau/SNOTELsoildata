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
octAirTmean = climData(:, 14);
octAirTsd = climData(:, 15);
novAirTmean = climData(:, 16);
novAirTsd = climData(:, 17);
decAirTmean = climData(:, 18);
decAirTsd = climData(:, 19);
janAirTmean = climData(:, 20);
janAirTsd = climData(:, 21);
febAirTmean = climData(:, 22);
febAirTsd = climData(:, 23);
marAirTmean = climData(:, 24);
marAirTsd = climData(:, 25);
aprAirTmean = climData(:, 26);
aprAirTsd = climData(:, 27);
mayAirTmean = climData(:, 28);
mayAirTsd = climData(:, 29);
junAirTmean = climData(:, 30);
junAirTsd = climData(:, 31);
julAirTmean = climData(:, 32);
julAirTsd = climData(:, 33);
augAirTmean = climData(:, 34);
augAirTsd = climData(:, 35);
sepAirTmean = climData(:, 36);
sepAirTsd = climData(:, 37);
maat = climData(:, 38);
sdAnnAirT = climData(:, 39);

preonsetAirT = climData(:, 40);
preonsetAirTsd = climData(:, 41);
premeltAirT = climData(:, 42);
premeltAirTsd = climData(:, 43);
postmeltAirT = climData(:, 44);
postmeltAirTsd = climData(:, 45);
elev = climData(:, 46);
lat = climData(:, 47);
lon = climData(:, 48);
ltMeanSWE = climData(:, 49);
ltMeanPrecip = climData(:, 50);

% Parse soilwatersummary
site_sw = vwcData(:, 1);
year_sw = vwcData(:, 2);
oct5cmSMmean = vwcData(:, 3);
oct5cmSMsd = vwcData(:, 4);
oct20cmSMmean = vwcData(:, 5);
oct20cmSMsd = vwcData(:, 6);
oct50cmSMmean = vwcData(:, 7);
oct50cmSMsd = vwcData(:, 8);

dec5cmSMmean = vwcData(:, 15);
dec5cmSMsd = vwcData(:, 16);
dec20cmSMmean = vwcData(:, 17);
dec20cmSMsd = vwcData(:, 18);
dec50cmSMmean = vwcData(:, 19);
dec50cmSMsd = vwcData(:, 20);

feb5cmSMmean = vwcData(:, 27);
feb5cmSMsd = vwcData(:, 28);
feb20cmSMmean = vwcData(:, 29);
feb20cmSMsd = vwcData(:, 30);
feb50cmSMmean = vwcData(:, 31);
feb50cmSMsd = vwcData(:, 32);

apr5cmSMmean = vwcData(:, 45);
apr5cmSMsd = vwcData(:, 46);
apr20cmSMmean = vwcData(:, 47);
apr20cmSMsd = vwcData(:, 48);
apr50cmSMmean = vwcData(:, 49);
apr50cmSMsd = vwcData(:, 50);

% These repeat through sept (end of wy)
ond5cmSMmean = vwcData(:, 73);
ond5cmSMsd = vwcData(:, 74);
ond20cmSMmean = vwcData(:, 75);
ond20cmSMsd = vwcData(:, 76);
ond50cmSMmean = vwcData(:, 77);
ond50cmSMsd = vwcData(:, 78);
jfm5cmSMmean = vwcData(:, 79);
jfm5cmSMsd = vwcData(:, 80);
jfm20cmSMmean = vwcData(:, 81);
jfm20cmSMsd = vwcData(:, 82);
jfm50cmSMmean = vwcData(:, 83);
jfm50cmSMsd = vwcData(:, 84);
amj5cmSMmean = vwcData(:, 85);
amj5cmSMsd = vwcData(:, 86);
amj20cmSMmean = vwcData(:, 87);
amj20cmSMsd = vwcData(:, 88);
amj50cmSMmean = vwcData(:, 89);
amj50cmSMsd = vwcData(:, 90);
jas5cmSMmean = vwcData(:, 91);
jas5cmSMsd = vwcData(:, 92);
jas20cmSMmean = vwcData(:, 93);
jas20cmSMsd = vwcData(:, 94);
jas50cmSMmean = vwcData(:, 95);
jas50cmSMsd = vwcData(:, 96);

preonset5cmSM = vwcData(:, 97);
preonset5cmSMsd = vwcData(:, 98);
preonset20cmSM = vwcData(:, 99);
preonset20cmSMsd = vwcData(:, 100);
preonset50cmSM = vwcData(:, 101);
preonset50cmSMsd = vwcData(:, 102);

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
ond5cmSTmean = tsData(:, 73);
ond5cmSTsd = tsData(:, 74);
ond20cmSTmean = tsData(:, 75);
ond20cmSTsd = tsData(:, 76);
ond50cmSTmean = tsData(:, 77);
ond50cmSTsd = tsData(:, 78);
jfm5cmSTmean = tsData(:, 79);
jfm5cmSTsd = tsData(:, 80);
jfm20cmSTmean = tsData(:, 81);
jfm20cmSTsd = tsData(:, 82);
jfm50cmSTmean = tsData(:, 83);
jfm50cmSTsd = tsData(:, 84);
amj5cmSTmean = tsData(:, 85);
amj5cmSTsd = tsData(:, 86);
amj20cmSTmean = tsData(:, 87);
amj20cmSTsd = tsData(:, 88);
amj50cmSTmean = tsData(:, 89);
amj50cmSTsd = tsData(:, 90);
jas5cmSTmean = tsData(:, 91);
jas5cmSTsd = tsData(:, 92);
jas20cmSTmean = tsData(:, 93);
jas20cmSTsd = tsData(:, 94);
jas50cmSTmean = tsData(:, 95);
jas50cmSTsd = tsData(:, 96);

% Seasonal/yearly soil temp means
mast5cm = tsData(:, 97);
sdast5cm = tsData(:, 98);
mast20cm = tsData(:, 99);
sdast20cm = tsData(:, 100);
mast50cm = tsData(:, 101);
sdast50cm = tsData(:, 102);

% Snowcovered soil temp means
snowcovMeanST5cm = tsData(:, 103);
snowcovStdST5cm = tsData(:, 104);
snowcovMeanST20cm = tsData(:, 105);
snowcovStdST20cm = tsData(:, 106);
snowcovMeanST50cm = tsData(:, 107);
snowcovStdST50cm = tsData(:, 108);

% Snowfree soil temp means
snowfreeMeanST5cm = tsData(:, 109);
snowfreeStdST5cm = tsData(:, 110);
snowfreeMeanST20cm = tsData(:, 111);
snowfreeStdST20cm = tsData(:, 112);
snowfreeMeanST50cm = tsData(:, 113);
snowfreeStdST50cm = tsData(:, 114);

preonset5cmST = tsData(:, 115);
preonset5cmSTsd = tsData(:, 116);
preonset20cmST = tsData(:, 117);
preonset20cmSTsd = tsData(:, 118);
preonset50cmST = tsData(:, 119);
preonset50cmSTsd = tsData(:, 120);
premelt5cmST = tsData(:, 121);
premelt5cmSTsd = tsData(:, 122);
postmelt5cmST = tsData(:, 123);
postmelt5cmSTsd = tsData(:, 124);

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
title('Day of snowmelt');


% %------------------------------------------------------------------
% % FIG 3 - 3d histogram of some data above (just as an example)
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
set(h, 'Name', 'MAST - MAT for all wateryears');
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

subplot (2, 2, 1)
% Set binning parameters
topEdge = 45; % define limits
botEdge = 0; % define limits
numBins = 10; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1));

x = preonset5cmSM; %split into x and y
y = oct5cmSMmean;
y2 = oct5cmSMsd;
%[binMean1, binMean2] = bin(x, y, y2);
[binMeanDec, binSdDec] = binseries(x, y, y2, topEdge, botEdge, numBins);
errorbar(xax(1:numBins), binMeanDec, binSdDec, '.g');
hold on;
x = preonset20cmSM; %split into x and y
y = oct20cmSMmean;
y2 = oct20cmSMsd;
%[binMean1, binMean2] = bin(x, y, y2);
[binMeanDec, binSdDec] = binseries(x, y, y2, topEdge, botEdge, numBins);
errorbar(xax(1:numBins), binMeanDec, binSdDec, '.b');
x = preonset50cmSM; %split into x and y
y = oct50cmSMmean;
y2 = oct50cmSMsd;
%[binMean1, binMean2] = bin(x, y, y2);
[binMeanDec, binSdDec] = binseries(x, y, y2, topEdge, botEdge, numBins);
errorbar(xax(1:numBins), binMeanDec, binSdDec, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Wateryear pre-snowpack soilVWC (%)');
ylabel('Mean VWC (%)');
xlim([-5, 45]);
ylim([-5, 45]);
title('Oct VWC vs pre-snowpack VWC');

subplot (2, 2, 2)
% Set binning parameters
topEdge = 45; % define limits
botEdge = 0; % define limits
numBins = 10; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1));

x = preonset5cmSM; %split into x and y
y = dec5cmSMmean;
y2 = dec5cmSMsd;
%[binMean1, binMean2] = bin(x, y, y2);
[binMeanDec, binSdDec] = binseries(x, y, y2, topEdge, botEdge, numBins);
errorbar(xax(1:numBins), binMeanDec, binSdDec, '.g');
hold on;
x = preonset20cmSM; %split into x and y
y = dec20cmSMmean;
y2 = dec20cmSMsd;
%[binMean1, binMean2] = bin(x, y, y2);
[binMeanDec, binSdDec] = binseries(x, y, y2, topEdge, botEdge, numBins);
errorbar(xax(1:numBins), binMeanDec, binSdDec, '.b');
x = preonset50cmSM; %split into x and y
y = dec50cmSMmean;
y2 = dec50cmSMsd;
%[binMean1, binMean2] = bin(x, y, y2);
[binMeanDec, binSdDec] = binseries(x, y, y2, topEdge, botEdge, numBins);
errorbar(xax(1:numBins), binMeanDec, binSdDec, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Wateryear pre-snowpack soilVWC (%)');
ylabel('Mean VWC (%)');
xlim([-5, 45]);
ylim([-5, 45]);
title('Dec VWC vs pre-snowpack VWC');

subplot (2, 2, 3)
% Set binning parameters
topEdge = 45; % define limits
botEdge = 0; % define limits
numBins = 10; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1));

x = preonset5cmSM; %split into x and y
y = feb5cmSMmean;
y2 = feb5cmSMsd;
%[binMean1, binMean2] = bin(x, y, y2);
[binMeanFeb,binSdFeb] = binseries(x, y, y2, topEdge, botEdge, numBins);
errorbar(xax(1:numBins), binMeanFeb, binSdFeb, '.g');
hold on;
x = preonset20cmSM; %split into x and y
y = feb20cmSMmean;
y2 = feb20cmSMsd;
%[binMean1, binMean2] = bin(x, y, y2);
[binMeanFeb,binSdFeb] = binseries(x, y, y2, topEdge, botEdge, numBins);
errorbar(xax(1:numBins), binMeanFeb, binSdFeb, '.b');
x = preonset50cmSM; %split into x and y
y = feb50cmSMmean;
y2 = feb50cmSMsd;
%[binMean1, binMean2] = bin(x, y, y2);
[binMeanFeb,binSdFeb] = binseries(x, y, y2, topEdge, botEdge, numBins);
errorbar(xax(1:numBins), binMeanFeb, binSdFeb, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Wateryear pre-snowpack soilVWC (%)');
ylabel('Mean VWC (%)');
xlim([-5, 45]);
ylim([-5, 45]);
title('February VWC vs pre-snowpack VWC');

subplot (2, 2, 4)
% Set binning parameters
topEdge = 45; % define limits
botEdge = 0; % define limits
numBins = 10; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1));

x = preonset5cmSM; %split into x and y
y = apr5cmSMmean;
y2 = apr5cmSMsd;
%[binMean1, binMean2] = bin(x, y, y2);
[binMeanApr,binSdApr] = binseries(x, y, y2, topEdge, botEdge, numBins);
errorbar(xax(1:numBins), binMeanApr, binSdApr, '.g');
hold on;
x = preonset20cmSM; %split into x and y
y = apr20cmSMmean;
y2 = apr20cmSMsd;
%[binMean1, binMean2] = bin(x, y, y2);
[binMeanApr,binSdApr] = binseries(x, y, y2, topEdge, botEdge, numBins);
errorbar(xax(1:numBins), binMeanApr, binSdApr, '.b');
x = preonset50cmSM; %split into x and y
y = apr50cmSMmean;
y2 = apr50cmSMsd;
%[binMean1, binMean2] = bin(x, y, y2);
[binMeanApr,binSdApr] = binseries(x, y, y2, topEdge, botEdge, numBins);
errorbar(xax(1:numBins), binMeanApr, binSdApr, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Wateryear pre-snowpack soilVWC (%)');
ylabel('Mean VWC (%)');
xlim([-5, 45]);
ylim([-5, 45]);
title('April VWC vs pre-snowpack VWC');

junk = 99;
end



