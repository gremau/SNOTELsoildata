% plot_snotel_summary.m
%
% Reads the outputs from summarize_wateryear.m and makes a number of plots
% characterizing variability in soil moisture and soil temperature in the
% SNOTEL network.
%
% Feb 20, 2012 - Greg Maurer

clear;          % clear memory
close all;      % clear any figures
fignum=0;       % used to increment figure number for plots

addpath('/home/greg/data/code_resources/m_common/'); % access to hist2, etc

% Set data path and file name, read in file
rawdatapath = '../rawdata/soilsensors_hourly/';
processeddatapath = '../processed_data/';

climatedata = csvread([processeddatapath 'wyear_climatesummary.csv']);
soilwaterdata = csvread([processeddatapath 'wyear_soilwatersummary.csv']);
normsoilwaterdata = csvread([processeddatapath ...
    'wyear_soilwatersummary_smnorm.csv']);
soiltempdata = csvread([processeddatapath 'wyear_soiltempsummary.csv']);
soilsiteyears = dlmread([rawdatapath 'sitelist.txt'], ',');
soilsites = unique(soilsiteyears(:, 1));

site_cl = climatedata(:, 1);
year_cl = climatedata(:, 2);
maxswe = climatedata(:, 3)*25.4;
maxsweday = climatedata(:, 4);
maxdepth = climatedata(:, 5);
onsetdatenum = climatedata(:, 6);
meltdatenum = climatedata(:, 7);
onsetdoy = climatedata(:, 8);
meltdoy = climatedata(:, 9);
snowduration = climatedata(:, 10);
accumprecip = climatedata(:, 11);
JASprecip = climatedata(:, 12);
octAirTmean = climatedata(:, 13);
octAirTsd = climatedata(:, 14);
novAirTmean = climatedata(:, 15);
novAirTsd = climatedata(:, 16);
decAirTmean = climatedata(:, 17);
decAirTsd = climatedata(:, 18);
janAirTmean = climatedata(:, 19);
janAirTsd = climatedata(:, 20);
febAirTmean = climatedata(:, 21);
febAirTsd = climatedata(:, 22);
marAirTmean = climatedata(:, 23);
marAirTsd = climatedata(:, 24);
aprAirTmean = climatedata(:, 25);
aprAirTsd = climatedata(:, 26);
mayAirTmean = climatedata(:, 27);
mayAirTsd = climatedata(:, 28);
junAirTmean = climatedata(:, 29);
junAirTsd = climatedata(:, 30);
julAirTmean = climatedata(:, 31);
julAirTsd = climatedata(:, 32);
augAirTmean = climatedata(:, 33);
augAirTsd = climatedata(:, 34);
sepAirTmean = climatedata(:, 35);
sepAirTsd = climatedata(:, 36);
meanAnnAirT = climatedata(:, 37);
elev = climatedata(:, 38);
lat = climatedata(:, 39);
lon = climatedata(:, 40);
ltMeanSWE = climatedata(:, 41);
ltMeanPrecip = climatedata(:, 42);

% Parse soilwatersummary
site_sw = soilwaterdata(:, 1);
year_sw = soilwaterdata(:, 2);
oct5cmSMmean = soilwaterdata(:, 3);
oct5cmSMsd = soilwaterdata(:, 4);
oct20cmSMmean = soilwaterdata(:, 5);
oct20cmSMsd = soilwaterdata(:, 6);
oct50cmSMmean = soilwaterdata(:, 7);
oct50cmSMsd = soilwaterdata(:, 8);

feb5cmSMmean = soilwaterdata(:, 27);
feb5cmSMsd = soilwaterdata(:, 28);
feb20cmSMmean = soilwaterdata(:, 29);
feb20cmSMsd = soilwaterdata(:, 30);
feb50cmSMmean = soilwaterdata(:, 31);
feb50cmSMsd = soilwaterdata(:, 32);

apr5cmSMmean = soilwaterdata(:, 45);
apr5cmSMsd = soilwaterdata(:, 46);
apr20cmSMmean = soilwaterdata(:, 47);
apr20cmSMsd = soilwaterdata(:, 48);
apr50cmSMmean = soilwaterdata(:, 49);
apr50cmSMsd = soilwaterdata(:, 50);


% These repeat through sept (end of wy)
ond5cmSMmean = soilwaterdata(:, 73);
ond5cmSMsd = soilwaterdata(:, 74);
ond20cmSMmean = soilwaterdata(:, 75);
ond20cmSMsd = soilwaterdata(:, 76);
ond50cmSMmean = soilwaterdata(:, 77);
ond50cmSMsd = soilwaterdata(:, 78);
jfm5cmSMmean = soilwaterdata(:, 79);
jfm5cmSMsd = soilwaterdata(:, 80);
jfm20cmSMmean = soilwaterdata(:, 81);
jfm20cmSMsd = soilwaterdata(:, 82);
jfm50cmSMmean = soilwaterdata(:, 83);
jfm50cmSMsd = soilwaterdata(:, 84);
amj5cmSMmean = soilwaterdata(:, 85);
amj5cmSMsd = soilwaterdata(:, 86);
amj20cmSMmean = soilwaterdata(:, 87);
amj20cmSMsd = soilwaterdata(:, 88);
amj50cmSMmean = soilwaterdata(:, 89);
amj50cmSMsd = soilwaterdata(:, 90);
jas5cmSMmean = soilwaterdata(:, 91);
jas5cmSMsd = soilwaterdata(:, 92);
jas20cmSMmean = soilwaterdata(:, 93);
jas20cmSMsd = soilwaterdata(:, 94);
jas50cmSMmean = soilwaterdata(:, 95);
jas50cmSMsd = soilwaterdata(:, 96);

preonsetAirT = soilwaterdata(:, 97);
preonset5cmSM = soilwaterdata(:, 98);
preonset20cmSM = soilwaterdata(:, 99);
preonset50cmSM = soilwaterdata(:, 100);
premeltAirT = soilwaterdata(:, 101);
postmeltAirT = soilwaterdata(:, 102);


% Parse NORMALIZED soilwatersummary
normamj5cmSMmean = normsoilwaterdata(:, 85);
normamj5cmSMsd = normsoilwaterdata(:, 86);
normamj20cmSMmean = normsoilwaterdata(:, 87);
normamj20cmSMsd = normsoilwaterdata(:, 88);
normamj50cmSMmean = normsoilwaterdata(:, 89);
normamj50cmSMsd = normsoilwaterdata(:, 90);
normjas5cmSMmean = normsoilwaterdata(:, 91);
normjas5cmSMsd = normsoilwaterdata(:, 92);
normjas20cmSMmean = normsoilwaterdata(:, 93);
normjas20cmSMsd = normsoilwaterdata(:, 94);
normjas50cmSMmean = normsoilwaterdata(:, 95);
normjas50cmSMsd = normsoilwaterdata(:, 96);


% Parse soiltempsummary
site_st = soiltempdata(:, 1);
year_st = soiltempdata(:, 2);
oct5cmSTmean = soiltempdata(:, 3);
oct5cmSTsd = soiltempdata(:, 4);
oct20cmSTmean = soiltempdata(:, 5);
oct20cmSTsd = soiltempdata(:, 6);
oct50cmSTmean = soiltempdata(:, 7);
oct50cmSTsd = soiltempdata(:, 8);
% These repeat through sept (end of wy)
ond5cmSTmean = soiltempdata(:, 73);
ond5cmSTsd = soiltempdata(:, 74);
ond20cmSTmean = soiltempdata(:, 75);
ond20cmSTsd = soiltempdata(:, 76);
ond50cmSTmean = soiltempdata(:, 77);
ond50cmSTsd = soiltempdata(:, 78);
jfm5cmSTmean = soiltempdata(:, 79);
jfm5cmSTsd = soiltempdata(:, 80);
jfm20cmSTmean = soiltempdata(:, 81);
jfm20cmSTsd = soiltempdata(:, 82);
jfm50cmSTmean = soiltempdata(:, 83);
jfm50cmSTsd = soiltempdata(:, 84);
amj5cmSTmean = soiltempdata(:, 85);
amj5cmSTsd = soiltempdata(:, 86);
amj20cmSTmean = soiltempdata(:, 87);
amj20cmSTsd = soiltempdata(:, 88);
amj50cmSTmean = soiltempdata(:, 89);
amj50cmSTsd = soiltempdata(:, 90);
jas5cmSTmean = soiltempdata(:, 91);
jas5cmSTsd = soiltempdata(:, 92);
jas20cmSTmean = soiltempdata(:, 93);
jas20cmSTsd = soiltempdata(:, 94);
jas50cmSTmean = soiltempdata(:, 95);
jas50cmSTsd = soiltempdata(:, 96);

% Seasonal/yearly soil temp means
mast5cm = soiltempdata(:, 97);
sdast5cm = soiltempdata(:, 98);
mast20cm = soiltempdata(:, 99);
sdast20cm = soiltempdata(:, 100);
mast50cm = soiltempdata(:, 101);
sdast50cm = soiltempdata(:, 102);

% Snowcovered soil temp means
snowcovMeanST5cm = soiltempdata(:, 103);
snowcovStdST5cm = soiltempdata(:, 104);
snowcovMeanST20cm = soiltempdata(:, 105);
snowcovStdST20cm = soiltempdata(:, 106);
snowcovMeanST50cm = soiltempdata(:, 107);
snowcovStdST50cm = soiltempdata(:, 108);

% Snowfree soil temp means
snowfreeMeanST5cm = soiltempdata(:, 109);
snowfreeStdST5cm = soiltempdata(:, 110);
snowfreeMeanST20cm = soiltempdata(:, 111);
snowfreeStdST20cm = soiltempdata(:, 112);
snowfreeMeanST50cm = soiltempdata(:, 113);
snowfreeStdST50cm = soiltempdata(:, 114);

preonset5cmST = soiltempdata(:, 115);
preonset20cmST = soiltempdata(:, 116);
preonset50cmST = soiltempdata(:, 117);
premelt5cmST = soiltempdata(:, 118);
postmelt5cmST = soiltempdata(:, 119);

% TEMPORARY outlier repair
test = mast5cm>15 | mast20cm>15 | mast50cm>15 | mast5cm<-3.5;
mast5cm(test) = nan; mast20cm(test) = nan; mast50cm(test) = nan;
test = meanAnnAirT>13 | meanAnnAirT<-2; 
meanAnnAirT(test) = nan;

% Get a subset of climatedata that corresponds with available soildata
matchtest = ismember(climatedata(:, 1:2), soilsiteyears(:, 1:2), 'rows');
matchsets = climatedata(matchtest, :);

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
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Wateryear snowpack metrics 1 - all sites');

colormap(jet);

subplot (2, 2, 1)
% 2D histogram
test = ~isnan(elev) & ~isnan(maxswe);
x = elev(test);
y = maxswe(test);
xedges = linspace(500, 4000, 60);
yedges = linspace(0, 2500, 60);
histmat = hist2(x, y, xedges, yedges);  % hist2 is from the matlab user forum
pcolor(xedges, yedges, histmat');
title('Peak SWE across the network');

subplot (2, 2, 2)
% 2D histogram
test = ~isnan(elev) & ~isnan(accumprecip);
x = elev(test);
y = accumprecip(test);
xedges = linspace(500, 4000, 60);
yedges = linspace(0, 120, 60);
histmat = hist2(x, y, xedges, yedges);  % hist2 is from the matlab user forum
pcolor(xedges, yedges, histmat');
title('Wateryear precip across the network');


subplot (2, 2, 3)
plot(elev, maxswe, 'ok', 'MarkerFaceColor', 'k');
hold on;
testsoil = ismember(site_cl, soilsites);
plot(elev(testsoil), maxswe(testsoil), 'ob');
xlabel('Elevation (m)');
ylabel('Peak SWE (mm)');
legend('Intermountain west', 'Soil sites');


subplot (2, 2, 4)
plot(elev, accumprecip, 'ob', 'MarkerFaceColor', 'b');
hold on;
plot(elev(testsoil), accumprecip(testsoil), 'ok');
xlabel('Elevation (m)');
ylabel('Wateryear precip (mm)');
legend('Intermountain west', 'Soil sites');


%  --------------------------------------------------------
% FIG 2 - Add MAT to analysis above
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Wateryear snowpack metrics 1 - all sites');

colormap(jet);

subplot (2, 1, 1)
% 2D histogram
test = ~isnan(elev) & ~isnan(meanAnnAirT);
x = elev(test);
y = meanAnnAirT(test);
xedges = linspace(500, 4000, 60);
yedges = linspace(-2, 14, 60);
histmat = hist2(x, y, xedges, yedges);  % hist2 is from the matlab user forum
pcolor(xedges, yedges, histmat');
title('Mean wateryear airT across network');

subplot (2, 1, 2)
plot(elev, meanAnnAirT, 'om', 'MarkerFaceColor', 'm');
hold on;
testsoil = ismember(site_cl, soilsites);
plot(elev(testsoil), meanAnnAirT(testsoil), 'ok');
xlabel('Elevation (m)');
ylabel('Mean wateryear airT (^oC)');
legend('Intermountain west', 'Soil sites');

%-------------------------------------------------------------
% FIG 3 
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Wateryear snowpack metrics 2 - all sites');

subplot (2, 2, 1)
plot(elev, onsetdoy, 'ob');
ylim([0 250]);
xlabel('Elevation');
ylabel('Day of water year');
title('Snowpack onset day');

subplot (2, 2, 2)
% 2D histogram
test = ~isnan(elev) & ~isnan(onsetdoy);
x = elev(test);
y = onsetdoy(test);
xedges = linspace(0, 4000, 100);
yedges = linspace(0, 250, 100);
histmat = hist2(x, y, xedges, yedges);  % hist2 is from the matlab user forum
pcolor(xedges, yedges, histmat'); 

subplot (2, 2, 3)
plot(elev, meltdoy, 'ob');
xlabel('Elevation');
ylabel('Day of water year');
title('Day of snowmelt');

subplot (2, 2, 4)
% 2D histogram
test = ~isnan(elev) & ~isnan(meltdoy);
x = elev(test);
y = meltdoy(test);
xedges = linspace(0, 4000, 100);
yedges = linspace(0, 400, 100);
histmat = hist2(x, y, xedges, yedges);  % hist2 is from the matlab user forum
pcolor(xedges, yedges, histmat');

%----------------------------------------------------
% FIG 4 - Wateryear precip metrics
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Wateryear precip metrics 1 - all sites');


% First four subplots on the left are for quarters
subplot (2, 2, 1)
plot(elev, accumprecip, 'ob');
hold on
%plot(elev, ltMeanSWE, 'or');
plot(elev(testsoil), accumprecip(testsoil), 'ok');
xlabel('Elevation');
ylabel('Precip (mm)');
legend('Individual years', 'Soil sites');
title('Total wateryear precip');

subplot (2, 2, 2)
% 2D histogram
test = ~isnan(elev) & ~isnan(accumprecip);
x = elev(test);
y = accumprecip(test);
xedges = linspace(0, 4000, 100);
yedges = linspace(0, 120, 100);
histmat = hist2(x, y, xedges, yedges);  % hist2 is from the matlab user forum
pcolor(xedges, yedges, histmat');

% First four subplots on the left are for quarters
subplot (2, 2, 3)
plot(elev, JASprecip, 'ob');
hold on
%plot(elev, ltMeanSWE, 'or');
plot(elev(testsoil), JASprecip(testsoil), 'ok');
xlabel('Elevation');
ylabel('Precip (mm)');
ylim([0, 120]);
legend('Individual years', 'Soil sites');
title('Summer (JAS) precip');

subplot (2, 2, 4)
% 2D histogram
test = ~isnan(elev) & ~isnan(JASprecip);
x = elev(test);
y = JASprecip(test);
xedges = linspace(0, 4000, 100);
yedges = linspace(0, 120, 100);
histmat = hist2(x, y, xedges, yedges);  % hist2 is from the matlab user forum
pcolor(xedges, yedges, histmat');

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

%----------------------------------------------------
% FIG 6 - Plot NORMALIZED soil moisture vs max swe  and snowmelt day
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Wateryear NORMALIZED seasonal soil moisture');

subplot (2, 2, 1)
plot(maxswe(matchtest), normamj5cmSMmean,  '.g', ...
    maxswe(matchtest), normamj20cmSMmean, '.b', ...
    maxswe(matchtest), normamj50cmSMmean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Peak SWE (mm)');
ylabel('VWC (%)');
title(' Apr, May, Jun VWC');

subplot (2, 2, 2)
plot(maxswe(matchtest), normjas5cmSMmean,  '.g', ...
    maxswe(matchtest), normjas20cmSMmean, '.b', ...
    maxswe(matchtest), normjas50cmSMmean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Peak SWE (mm)');
ylabel('VWC (%)');
title(' July, Aug, Sept VWC');

subplot (2, 2, 3)
plot(meltdoy(matchtest), normamj5cmSMmean,  '.g', ...
    meltdoy(matchtest), normamj20cmSMmean, '.b', ...
    meltdoy(matchtest), normamj50cmSMmean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Day of wateryear');
ylabel('VWC (%)');
title(' Apr, May, Jun VWC vs snowmelt day');

subplot (2, 2, 4)
plot(meltdoy(matchtest), normjas5cmSMmean,  '.g', ...
    meltdoy(matchtest), normjas20cmSMmean, '.b', ...
    meltdoy(matchtest), normjas50cmSMmean, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Day of wateryear (mm)');
ylabel('VWC (%)');
title(' July, Aug, Sept VWC vs snowmeltday');

%--------------------------------------------------------------
% FIG 7 - Plot MAST % MAT vs max swe, snowcover duration, elevation
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Compare MAST & MAT (vs maxSWE, snowduration, elev');
subplot (2, 2, 1)
plot(maxswe(matchtest), meanAnnAirT(matchtest), 'ok', ...
    'MarkerFaceColor', [0.7 0.7 0.7]);
hold on;
plot(maxswe(matchtest), mast5cm, '.b', ...
    maxswe(matchtest), mast20cm, '.b', ...
    maxswe(matchtest), mast50cm, '.b');
legend('Mean Air T', 'Mean Soil T');
xlabel('Peak SWE (mm)');
ylabel('^oC');
title('Mean wateryear AirT & SoilT vs peak SWE');

subplot (2, 2, 2)
plot(elev(matchtest), meanAnnAirT(matchtest), 'ok', ...
    'MarkerFaceColor', [0.7 0.7 0.7]);
hold on;
plot(elev(matchtest), mast5cm, '.b', ...
    elev(matchtest), mast20cm, '.b',...
    elev(matchtest), mast50cm, '.b');
legend('Mean Air T', 'Mean Soil T');
xlabel('Elevation (m)');
ylabel('^oC');
title('Mean wateryear AirT & SoilT vs Elevation');

subplot (2, 2, 3)
plot(snowduration(matchtest), meanAnnAirT(matchtest), 'ok', ...
    'MarkerFaceColor', [0.7 0.7 0.7]);
hold on;
plot(snowduration(matchtest), mast5cm, '.b', ...
    snowduration(matchtest), mast20cm, '.b',...
    snowduration(matchtest), mast50cm, '.b');
legend('Mean Air T', 'Mean Soil T');
xlabel('Snowpack duration (days)');
ylabel('^oC');
title('... vs Snowpack duration');

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
y = meanAnnAirT(matchtest);
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
y = meanAnnAirT(matchtest);
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
y = meanAnnAirT(matchtest);
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
plot(meanAnnAirT(matchtest), mast5cm, '.g', ...
    meanAnnAirT(matchtest), mast20cm, '.b',...
    meanAnnAirT(matchtest), mast50cm, '.k');
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
matchMAT = meanAnnAirT(matchtest);

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
plot(preonsetAirT, snowcovMeanST5cm, '.g', ...
    preonsetAirT, snowcovMeanST20cm, '.b', ...
    preonsetAirT, snowcovMeanST50cm, '.k');
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
plot(preonset5cmSM, ond5cmSMmean, '.g', ...
    preonset20cmSM, ond20cmSMmean, '.b', ...
    preonset50cmSM, ond50cmSMmean, '.k');
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
numBins = 15; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1)) + 1.5;

x = preonset5cmSM; %split into x and y
y = ond5cmSMmean;
y2 = ond5cmSMsd;
%[binMean1, binMean2] = bin(x, y, y2);
[binMeanOND, binSdOND] = binseries(x, y, y2, topEdge, botEdge, numBins);
errorbar(xax(1:numBins), binMeanOND, binSdOND, '.g');
hold on;
x = preonset20cmSM; %split into x and y
y = ond20cmSMmean;
y2 = ond20cmSMsd;
%[binMean1, binMean2] = bin(x, y, y2);
[binMeanOND, binSdOND] = binseries(x, y, y2, topEdge, botEdge, numBins);
errorbar(xax(1:numBins), binMeanOND, binSdOND, '.b');
x = preonset50cmSM; %split into x and y
y = ond50cmSMmean;
y2 = ond50cmSMsd;
%[binMean1, binMean2] = bin(x, y, y2);
[binMeanOND, binSdOND] = binseries(x, y, y2, topEdge, botEdge, numBins);
errorbar(xax(1:numBins), binMeanOND, binSdOND, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Wateryear pre-snowpack soilVWC (%)');
ylabel('Mean VWC (%)');
title('Oct, Nov, Dec VWC vs pre-snowpack VWC');

subplot (2, 2, 2)
% Set binning parameters
topEdge = 45; % define limits
botEdge = 0; % define limits
numBins = 15; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1)) + 1.5;

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
title('February VWC vs pre-snowpack VWC');

subplot (2, 2, 3)
% Set binning parameters
topEdge = 45; % define limits
botEdge = 0; % define limits
numBins = 15; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1)) + 1.5;

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
ylim([0, 45]);
title('April VWC vs pre-snowpack VWC');



