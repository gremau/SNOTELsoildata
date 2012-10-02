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

climatedata = csvread([processeddatapath 'wyear_climatesummary.txt']);
soilwaterdata = csvread([processeddatapath 'wyear_soilwatersummary_hourly.txt']);
normsoilwaterdata = csvread([processeddatapath ...
    'wyear_soilwatersummary_hourly_smnorm.txt']);
soiltempdata = csvread([processeddatapath 'wyear_soiltempsummary.csv']);
soilsiteyears = dlmread([rawdatapath 'filelist.txt'], ',');
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

aug5cmSMmean = soilwaterdata(:, 63);
aug5cmSMsd = soilwaterdata(:, 64);
aug20cmSMmean = soilwaterdata(:, 65);
aug20cmSMsd = soilwaterdata(:, 66);
aug50cmSMmean = soilwaterdata(:, 67);
aug50cmSMsd = soilwaterdata(:, 68);

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

% preonsetAirT = soilwaterdata(:, 97);
% preonset5cmSM = soilwaterdata(:, 98);
% preonset20cmSM = soilwaterdata(:, 99);
% preonset50cmSM = soilwaterdata(:, 100);
% premeltAirT = soilwaterdata(:, 101);
% postmeltAirT = soilwaterdata(:, 102);


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


%----------------------------------------------------
% FIG 3 - Plot soil moisture vs max swe  and snowmelt day
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'JAS soil moisture vs snowmelt date');
% Set binning parameters
topEdge = 300; % define limits
botEdge = 90; % define limits
numBins = 15; % define number of bins
% And an xaxis to use
xax = (linspace(botEdge, topEdge, numBins+1))+5.25;

% 5cm
x = meltdoy(matchtest); %split into x and y
y = jas5cmSMmean;
y2 = jas5cmSMsd;
%[binMean1, binMean2] = bin(x, y, y2);
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins)

subplot(221);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([80 310]);
%ylim([0 1]);
xlabel('Snowmelt dowy');
title('5cm VWC');

% 20cm
y = jas20cmSMmean;
y2 = jas20cmSMsd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins)

subplot(222);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([80 310]);
%ylim([0 1]);
xlabel('Snowmelt dowy');
title('20cm VWC');

% 50cm
y = jas50cmSMmean;
y2 = jas50cmSMsd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins)

subplot(223);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([80 310]);
%ylim([0 1]);
xlabel('Snowmelt dowy');
title('50cm VWC');

%----------------------------------------------------
% FIG 4 - NORMALIZED JAS soil moisture vs snowmelt date
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'NORMALIZED JAS soil moisture vs snowmelt date');

% 5cm
x = meltdoy(matchtest); %split into x and y
y = normjas5cmSMmean;
y2 = normjas5cmSMsd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins)

subplot(221);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([80 310]);
ylim([0 1]);
xlabel('Snowmelt dowy');
title('5cm VWC');

% 20cm
y = normjas20cmSMmean;
y2 = normjas20cmSMsd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins)

subplot(222);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([80 310]);
ylim([0 1]);
xlabel('Snowmelt dowy');
title('20cm VWC');

% 50cm
y = normjas50cmSMmean;
y2 = normjas50cmSMsd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins)

subplot(223);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([80 310]);
ylim([0 1]);
xlabel('Snowmelt dowy');
title('50cm VWC');

%----------------------------------------------------
% FIG 5 - Plot NORMALIZED JAS soil moisture vs max swe
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'NORMALIZED JAS soil moisture vs max swe');

% First change some binning/plotting parameters
topEdge=1000;
xax = (linspace(botEdge, topEdge, numBins+1))+20;

% 5cm
x = maxswe(matchtest); %split into x and y
y = normjas5cmSMmean;
y2 = normjas5cmSMsd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins)

subplot(221);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([80 1000]);
ylim([0 1]);
xlabel('Peak SWE');
title('5cm VWC');

% 20cm
y = normjas20cmSMmean;
y2 = normjas20cmSMsd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins)

subplot(222);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([80 1000]);
ylim([0 1]);
xlabel('Peak SWE');
title('20cm VWC');

% 50cm
y = normjas50cmSMmean;
y2 = normjas50cmSMsd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins)

subplot(223);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([80 1000]);
ylim([0 1]);
xlabel('Peak SWE');
title('50cm VWC');

%----------------------------------------------------
% FIG 6 - Plot August soil moisture vs snowmelt date
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Aug soil moisture vs snomelt date');

% First change some binning/plotting parameters
topEdge=300;
xax = (linspace(botEdge, topEdge, numBins+1))+5.25;

% 5cm

x = meltdoy(matchtest); %split into x and y
y = aug5cmSMmean;
y2 = aug5cmSMsd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins)


subplot(221);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([80 310]);
%ylim([0 1]);
xlabel('Snowmelt dowy');
title('5cm VWC');

% 20cm
y = aug20cmSMmean;
y2 = aug20cmSMsd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins)

subplot(222);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([80 310]);
%ylim([0 1]);
xlabel('Snowmelt dowy');
title('20cm VWC');

% 50cm
y = aug50cmSMmean;
y2 = aug50cmSMsd;
[binMean, binStd] = binseries(x, y, y2, topEdge, botEdge, numBins)

subplot(223);
errorbar(xax(1:numBins), binMean, binStd, 'ok');
xlim([80 310]);
%ylim([0 1]);
xlabel('Snowmelt dowy');
title('50cm VWC');

end