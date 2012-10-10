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
addpath('/home/greg/data/code_resources/m_common/hline_vline/'); 

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

function [coefficients, rsq] = getstats(x, y);
    nantest = isnan(x) | isnan(y);
    x = x(~nantest);
    y = y(~nantest);
    coefficients = polyfit(x, y, 1);
    yfit = polyval(coefficients, x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
end

%-----PLOTTING----------------------------------------------------
% First define the 2d histogram plots
% colormap(jet)
% function @f = datadensity(xvar, yvar, xmax, ymax)
%     test = ~isnan(xvar) & ~isnan(yvar);
%     x = xvar(test);
%     y = yvar(test);
%     xedges = linspace(0, xmax, 75);
%     yedges = linspace(0, ymax, 75);
%     histmat = hist2(x, y, xedges, yedges);  % hist2 is from the matlab user forum
%     pcolor(xedges, yedges, histmat');
% end

%------------------------------------------------------------------
% FIG 1 - 
% fignum = fignum+1;
% h = figure(fignum);
% set(h, 'Name', 'Wateryear snowpack metrics 1 - all sites');
% 
% colormap(jet);
% 
% subplot (2, 2, 1)
% % 2D histogram
% test = ~isnan(elev) & ~isnan(maxswe);
% x = elev(test);
% y = maxswe(test);
% xedges = linspace(500, 4000, 60);
% yedges = linspace(0, 2500, 60);
% histmat = hist2(x, y, xedges, yedges);  % hist2 is from the matlab user forum
% pcolor(xedges, yedges, histmat');
% title('Peak SWE across the network');
% 
% subplot (2, 2, 2)
% % 2D histogram
% test = ~isnan(elev) & ~isnan(accumprecip);
% x = elev(test);
% y = accumprecip(test);
% xedges = linspace(500, 4000, 60);
% yedges = linspace(0, 120, 60);
% histmat = hist2(x, y, xedges, yedges);  % hist2 is from the matlab user forum
% pcolor(xedges, yedges, histmat');
% title('Wateryear precip across the network');
% 
% subplot (2, 2, 3)
% plot(elev, maxswe, 'ok', 'MarkerFaceColor', 'k');
% hold on;
% testsoil = ismember(site_cl, soilsites);
% plot(elev(testsoil), maxswe(testsoil), 'ob');
% xlabel('Elevation (m)');
% ylabel('Peak SWE (mm)');
% legend('Intermountain west', 'Soil sites');
% 
% subplot (2, 2, 4)
% plot(elev, accumprecip, 'ob', 'MarkerFaceColor', 'b');
% hold on;
% plot(elev(testsoil), accumprecip(testsoil), 'ok');
% xlabel('Elevation (m)');
% ylabel('Wateryear precip (mm)');
% legend('Intermountain west', 'Soil sites');


% FIG 2 - 
testsoil = ismember(site_cl, soilsites);
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Wateryear snowpack metrics 1 - all sites');

subplot (3, 2, 1)
xedges = linspace(500, 4000, 60);
networkhist = histc(elev, xedges);
soilhist = histc(elev(testsoil), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'r');
vline(nanmean(elev), ':k');
vline(nanmean(elev(testsoil)), ':r');
title('Elevation');

subplot (3, 2, 2)
xedges = linspace(-5, 20, 60);
networkhist = histc(maat, xedges);
soilhist = histc(maat(testsoil), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'r');
vline(nanmean(maat), ':k');
vline(nanmean(maat(testsoil)), ':r');
title('Mean wateryear air temp');

subplot (3, 2, 3)
xedges = linspace(0, 2500, 60);
networkhist = histc(accumprecip, xedges);
soilhist = histc(accumprecip(testsoil), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'r');
vline(nanmean(accumprecip), ':k');
vline(nanmean(accumprecip(testsoil)), ':r');
title('Wateryear precip');

subplot (3, 2, 4)
xedges = linspace(100, 2000, 60);
networkhist = histc(maxswe, xedges);
soilhist = histc(maxswe(testsoil), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'r');
vline(nanmean(maxswe), ':k');
vline(nanmean(maxswe(testsoil)), ':r');
title('Peak SWE');

subplot (3, 2, 5)
xedges = linspace(0, 365, 60);
networkhist = histc(totaldaysSC, xedges);
soilhist = histc(totaldaysSC(testsoil), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'r');
vline(nanmean(totaldaysSC), ':k');
vline(nanmean(totaldaysSC(testsoil)), ':r');
title('Total snowcovered days');

subplot (3, 2, 6)
xedges = linspace(0, 365, 60);
networkhist = histc(meltdoy, xedges);
soilhist = histc(meltdoy(testsoil), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'r');
vline(nanmean(meltdoy), ':k');
vline(nanmean(meltdoy(testsoil)), ':r');
title('Day of snowmelt');
% %  --------------------------------------------------------
% % FIG 2 - Add MAT to analysis above
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

% %-------------------------------------------------------------
% % FIG 3 
% fignum = fignum+1;
% h = figure(fignum);
% set(h, 'Name', 'Wateryear snowpack metrics 2 - all sites');
% 
% subplot (2, 2, 1)
% plot(elev, onsetdoy, 'ob');
% ylim([0 250]);
% xlabel('Elevation');
% ylabel('Day of water year');
% title('Snowpack onset day');
% 
% subplot (2, 2, 2)
% % 2D histogram
% test = ~isnan(elev) & ~isnan(onsetdoy);
% x = elev(test);
% y = onsetdoy(test);
% xedges = linspace(0, 4000, 100);
% yedges = linspace(0, 250, 100);
% histmat = hist2(x, y, xedges, yedges);  % hist2 is from the matlab user forum
% pcolor(xedges, yedges, histmat'); 
% 
% subplot (2, 2, 3)
% plot(elev, meltdoy, 'ob');
% xlabel('Elevation');
% ylabel('Day of water year');
% title('Day of snowmelt');
% 
% subplot (2, 2, 4)
% % 2D histogram
% test = ~isnan(elev) & ~isnan(meltdoy);
% x = elev(test);
% y = meltdoy(test);
% xedges = linspace(0, 4000, 100);
% yedges = linspace(0, 400, 100);
% histmat = hist2(x, y, xedges, yedges);  % hist2 is from the matlab user forum
% pcolor(xedges, yedges, histmat');

%----------------------------------------------------
% FIG 4 - Wateryear precip metrics
% fignum = fignum+1;
% h = figure(fignum);
% set(h, 'Name', 'Wateryear precip metrics 1 - all sites');
% 
% % First four subplots on the left are for quarters
% subplot (2, 2, 1)
% plot(elev, accumprecip, 'ob');
% hold on
% %plot(elev, ltMeanSWE, 'or');
% plot(elev(testsoil), accumprecip(testsoil), 'ok');
% xlabel('Elevation');
% ylabel('Precip (mm)');
% legend('Individual years', 'Soil sites');
% title('Total wateryear precip');
% 
% subplot (2, 2, 2)
% % 2D histogram
% test = ~isnan(elev) & ~isnan(accumprecip);
% x = elev(test);
% y = accumprecip(test);
% xedges = linspace(0, 4000, 100);
% yedges = linspace(0, 120, 100);
% histmat = hist2(x, y, xedges, yedges);  % hist2 is from the matlab user forum
% pcolor(xedges, yedges, histmat');
% 
% % First four subplots on the left are for quarters
% subplot (2, 2, 3)
% plot(elev, JASprecip, 'ob');
% hold on
% %plot(elev, ltMeanSWE, 'or');
% plot(elev(testsoil), JASprecip(testsoil), 'ok');
% xlabel('Elevation');
% ylabel('Precip (mm)');
% ylim([0, 120]);
% legend('Individual years', 'Soil sites');
% title('Summer (JAS) precip');
% 
% subplot (2, 2, 4)
% % 2D histogram
% test = ~isnan(elev) & ~isnan(JASprecip);
% x = elev(test);
% y = JASprecip(test);
% xedges = linspace(0, 4000, 100);
% yedges = linspace(0, 120, 100);
% histmat = hist2(x, y, xedges, yedges);  % hist2 is from the matlab user forum
% pcolor(xedges, yedges, histmat');

%----------------------------------------------------
% FIG 5 - Plot soil moisture vs max swe snowmelt day
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

%--------------------------------------------------------------
% FIG 7 - Plot MAST % MAT vs max swe, snowcover duration, elevation
% fignum = fignum+1;
% h = figure(fignum);
% set(h, 'Name', 'Compare MAST & MAT (vs maxSWE, snowduration, elev');
% subplot (2, 2, 1)
% plot(maxswe(matchtest), maat(matchtest), 'ok', ...
%     'MarkerFaceColor', [0.7 0.7 0.7]);
% hold on;
% plot(maxswe(matchtest), mast5cm, '.b', ...
%     maxswe(matchtest), mast20cm, '.b', ...
%     maxswe(matchtest), mast50cm, '.b');
% legend('Mean Air T', 'Mean Soil T');
% xlabel('Peak SWE (mm)');
% ylabel('^oC');
% title('Mean wateryear AirT & SoilT vs peak SWE');
% 
% subplot (2, 2, 2)
% plot(elev(matchtest), maat(matchtest), 'ok', ...
%     'MarkerFaceColor', [0.7 0.7 0.7]);
% hold on;
% plot(elev(matchtest), mast5cm, '.b', ...
%     elev(matchtest), mast20cm, '.b',...
%     elev(matchtest), mast50cm, '.b');
% legend('Mean Air T', 'Mean Soil T');
% xlabel('Elevation (m)');
% ylabel('^oC');
% title('Mean wateryear AirT & SoilT vs Elevation');
% 
% subplot (2, 2, 3)
% plot(snowduration(matchtest), maat(matchtest), 'ok', ...
%     'MarkerFaceColor', [0.7 0.7 0.7]);
% hold on;
% plot(snowduration(matchtest), mast5cm, '.b', ...
%     snowduration(matchtest), mast20cm, '.b',...
%     snowduration(matchtest), mast50cm, '.b');
% legend('Mean Air T', 'Mean Soil T');
% xlabel('Snowpack duration (days)');
% ylabel('^oC');
% title('... vs Snowpack duration');

%--------------------------------------------------------------
% FIG 7 - Plot MAST % MAT vs max swe, snowcover duration, elevation
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Compare MAST & MAT (vs maxSWE, snowduration, elev');
subplot (2, 2, 1)
plot(snowduration(matchtest), maat(matchtest), 'ok', ...
    'MarkerFaceColor', [0.7 0.7 0.7]);
hold on;
plot(snowduration(matchtest), mast5cm, '.b', ...
    snowduration(matchtest), mast20cm, '.b',...
    snowduration(matchtest), mast50cm, '.b');
legend('Mean Air T', 'Mean Soil T');
xlabel('Snowpack duration (days)');
ylabel('^oC');
title('MAST & MAT vs Snowpack duration');

subplot (2, 2, 2)
% plot(snowduration(matchtest), maat(matchtest), 'ok', ...
%     'MarkerFaceColor', [0.7 0.7 0.7]);
% hold on;
plot(snowduration(matchtest), mast5cm-maat(matchtest), '.r', ...
    snowduration(matchtest), mast20cm-maat(matchtest), '.g',...
    snowduration(matchtest), mast50cm-maat(matchtest), '.b');
legend('Mean Air T', 'Mean Soil T');
xlabel('Snowpack duration (days)');
ylabel('^oC');
title('... vs Snowpack duration');

bigtest = elev(matchtest)<2500 & elev(matchtest)>2000;
eltest = elev > 2000 & elev < 2500;
subplot (2, 2, 3)
% plot(snowduration(matchtest), maat(matchtest), 'ok', ...
%     'MarkerFaceColor', [0.7 0.7 0.7]);
% hold on;
plot(snowduration(matchtest & eltest), mast5cm(bigtest)-maat(matchtest & eltest), '.r', ...
    snowduration(matchtest & eltest), mast5cm(bigtest)-maat(matchtest & eltest), '.g',...
    snowduration(matchtest & eltest), mast5cm(bigtest)-maat(matchtest & eltest), '.b');
legend('Mean Air T', 'Mean Soil T');
xlabel('Snowpack duration (days)');
ylabel('^oC');
title('... vs Snowpack duration');

%--------------------------------------------------------------
% FIG  - Plot MAST vs snowpack duration for 3 sites
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Regress mast vs snowpack duration for 3 sites');
test = site_cl==828
site1dur = snowduration(matchtest & test);
site1maat = maat(matchtest & test);
site1mast = mast20cm(site_st==828);
test = site_cl==393
site2dur = snowduration(matchtest & test);
site2maat = maat(matchtest & test);
site2mast = mast20cm(site_st==393);
test = site_cl==582
site3dur = snowduration(matchtest & test);
site3mast = mast20cm(site_st==582);
site3maat = maat(matchtest & test);

subplot (3, 2, 1)
plot(site1dur, site1mast, '.b');
hold on;
[coefficients1, rsq1] = getstats(site1dur, site1mast);
linear_fit1 = polyval(coefficients1, [150 300]);
plot([150 300],linear_fit1,':k');
text(155, 4, ['r^2 = ' num2str(rsq1, 2)]); % r^2 values
xlabel('Snowpack duration (days)');
ylabel('MAST');
title('Trial Lake');

subplot (3, 2, 3)
plot(site2dur, site2mast, '.b');
hold on;
[coefficients2, rsq2] = getstats(site2dur, site2mast);
linear_fit2 = polyval(coefficients2, [150 300]);
plot([150 300],linear_fit2,':k');
text(155, 5, ['r^2 = ' num2str(rsq2, 2)]); % r^2 values
xlabel('Snowpack duration (days)');
ylabel('MAST');
title('Chalk Creek 1');

subplot (3, 2, 5)
plot(site3dur, site3mast, '.b');
hold on;
[coefficients3, rsq3] = getstats(site3dur, site3mast);
linear_fit3 = polyval(coefficients3, [150 300]);
plot([150 300],linear_fit3,':k');
text(155, 5, ['r^2 = ' num2str(rsq3, 2)]); % r^2 values
xlabel('Snowpack duration (days)');
ylabel('MAST');
title('Little Bear');

subplot (3, 2, 2)
plot(site1maat, site1mast, '.b')
hold on;
[coefficients4, rsq4] = getstats(site1maat, site1mast);
linear_fit4 = polyval(coefficients4, [0 9]);
plot([0 9],linear_fit4,':k');
text(3, 4, ['r^2 = ' num2str(rsq4, 2)]); % r^2 values
xlabel('MAT');
ylabel('MAST');
title('Trial Lake');

subplot (3, 2, 4)
plot(site2maat, site2mast, '.b')
hold on;
[coefficients5, rsq5] = getstats(site2maat, site2mast);
linear_fit5 = polyval(coefficients5, [0 9]);
plot([0 9],linear_fit4,':k');
text(6, 4, ['r^2 = ' num2str(rsq5, 2)]); % r^2 values
xlabel('MAT');
ylabel('MAST');
title('Chalk Creek 1');

subplot (3, 2, 6)
hold on;
[coefficients6, rsq6] = getstats(site3maat, site3mast);
linear_fit6 = polyval(coefficients6, [0 9]);
plot([0 9],linear_fit6,':k');
plot(site3maat, site3mast, '.b')
text(7, 7, ['r^2 = ' num2str(rsq6, 2)]); % r^2 values
xlabel('MAT');
ylabel('MAST');
title('Little Bear');

%--------------------------------------------------------------
% FIG 8 - Plot MAST % MAT vs max swe, snowcover duration, elevation
% SAME as above but binned using only 20cm data
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Compare MAST & MAT (vs maxSWE, snowduration, elev');

subplot (2, 2, 1)
% Set binning parameters
topEdge = 2000; % define limits
botEdge = 0; % define limits
numBins = 15; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1)) + 66;

x = maxswe(matchtest); %split into x and y
y = maat(matchtest);
y2 = mast20cm;
%[binMean1, binMean2] = bin(x, y, y2);
[binMeanAir, ~] = binseries(x, y, y2, topEdge, botEdge, numBins);
plot(xax(1:numBins), binMeanAir, 'ok', ...
    'MarkerFaceColor', [0.7 0.7 0.7]);
hold on;
y = mast20cm;
y2 = sdast20cm;
%[binMean1, binMean2] = bin(x, y, y2);
[binMeanSoil, binSdSoil] = binseries(x, y, y2, topEdge, botEdge, numBins);
errorbar(xax(1:numBins), binMeanSoil, binSdSoil, '.b');
legend('Mean Air T', 'Mean Soil T');
xlabel('Peak SWE (mm)');
ylabel('^oC');
title('Mean wateryear AirT & SoilT vs peak SWE');

subplot (2, 2, 2)
% Set binning parameters
topEdge = 3500; % define limits
botEdge = 920; % define limits
numBins = 15; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1)) + 85;

x = elev(matchtest);
y = maat(matchtest);
[binMeanAir, ~] = binseries(x, y, y2, topEdge, botEdge, numBins);
plot(xax(1:numBins), binMeanAir, 'ok', ...
    'MarkerFaceColor', [0.7 0.7 0.7]);
hold on;
y = mast20cm;
y2 = sdast20cm;
%[binMean1, binMean2] = bin(x, y, y2);
[binMeanSoil, binSdSoil] = binseries(x, y, y2, topEdge, botEdge, numBins);
errorbar(xax(1:numBins), binMeanSoil, binSdSoil, '.b');
legend('Mean Air T', 'Mean Soil T');
xlabel('Elevation (m)');
ylabel('^oC');
title('Mean wateryear AirT & SoilT vs Elevation');

subplot (2, 2, 3)
% Set binning parameters
topEdge = 300; % define limits
botEdge = 0; % define limits
numBins = 15; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1)) + 10;

x = snowduration(matchtest);
y = maat(matchtest);
[binMeanAir, ~] = binseries(x, y, y2, topEdge, botEdge, numBins);
plot(xax(1:numBins), binMeanAir, 'ok', ...
    'MarkerFaceColor', [0.7 0.7 0.7]);
hold on;
y = mast20cm;
y2 = sdast20cm;
%[binMean1, binMean2] = bin(x, y, y2);
[binMeanSoil, binSdSoil] = binseries(x, y, y2, topEdge, botEdge, numBins);
errorbar(xax(1:numBins), binMeanSoil, binSdSoil, '.b');
legend('Mean Air T', 'Mean Soil T');
xlabel('Snowpack duration (days)');
ylabel('^oC');
title('... vs Snowpack duration');

%--------------------------------------------------------------
% FIG 9 - Plot MAST (3 depths) vs SWE, snowcover duration, elevation, MAT
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'MAST vs SWE, snowduration, MAT, elevation');

subplot (2, 2, 1)
plot(maxswe(matchtest), mast5cm, '.g', ...
    maxswe(matchtest), mast20cm, '.b',...
    maxswe(matchtest), mast50cm, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Peak SWE (mm)');
ylabel('MWYST (^oC)');
title('Mean wateryear SoilT vs peak SWE');

subplot (2, 2, 2)
plot(maat(matchtest), mast5cm, '.g', ...
    maat(matchtest), mast20cm, '.b',...
    maat(matchtest), mast50cm, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Mean wyr AirT (^oC)');
ylabel('MWYST(^oC)');
title('Mean wateryear SoilT vs Mean wateryear AirT');

subplot (2, 2, 3)
plot(snowduration(matchtest), mast5cm, '.g', ...
    snowduration(matchtest), mast20cm, '.b',...
    snowduration(matchtest), mast50cm, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Snowpack duration (days)');
ylabel('MWYST (^oC)');
title('... vs snowpack duration');

subplot (2, 2, 4)
plot(elev(matchtest), mast5cm, '.g', ...
    elev(matchtest), mast20cm, '.b',...
    elev(matchtest), mast50cm, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Elevation (m)');
ylabel('MWYST (^oC)');
title('... vs Elevation');

%--------------------------------------------------------------
% FIG 10 - Plot MAST vs SWE, snowcover duration, MAT in elevation bins
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', ' MAST (50cm) vs SWE, MAT, snowcover duration, in elevation bins');

testhi = elev(matchtest) > 3000;
testmid = elev(matchtest) > 2500 & elev(matchtest) < 3000;
testlo = elev(matchtest) > 2000 & elev(matchtest) < 2500;

matchswe = maxswe(matchtest);
matchdur = snowduration(matchtest);
matchMAT = maat(matchtest);

subplot (2, 2, 1)
plot(matchswe(testhi), mast50cm(testhi), '.b', ...
    matchswe(testmid), mast50cm(testmid), '.g', ...
    matchswe(testlo), mast50cm(testlo), '.r');
legend('3000+', '2500-3000', '2000-2500cm');
xlabel('Peak SWE (mm)');
ylabel('MWYST (^oC)');
title('Mean wateryear SoilT vs peak SWE');

subplot (2, 2, 2)
plot(matchMAT(testhi), mast50cm(testhi), '.b', ...
    matchMAT(testmid), mast50cm(testmid), '.g',...
    matchMAT(testlo), mast50cm(testlo), '.r');
legend('3000+', '2500-3000', '2000-2500cm');
xlabel('Mean wateryear AirT (^oC)');
ylabel('MWYST (^oC)');
title('Mean wateryear SoilT vs. Mean wateryear AirT');

subplot (2, 2, 3)
plot(matchdur(testhi), mast50cm(testhi), '.b', ...
    matchdur(testmid), mast50cm(testmid), '.g',...
    matchdur(testlo), mast50cm(testlo), '.r');
legend('3000+', '2500-3000', '2000-2500cm');
xlabel('Snowpack duration (days)');
ylabel('MWYST(^oC)');
title('... vs Mean wateryear Air T');



%----------------------------------------------------
% FIG 11 - Plot snowcovered temp vs onset temps
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

%----------------------------------------------------
% FIG 12 - Plot winter soil moisture vs onset soil moisture
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
% FIG 13 - Plot winter soil moisture vs onset soil moisture
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



