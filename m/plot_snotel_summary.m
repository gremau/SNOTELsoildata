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
addpath('/home/greg/data/code_resources/m_common/nanstats/');
addpath('/home/greg/data/code_resources/m_common/');
addpath('/home/greg/data/code_resources/m_common/hline_vline/');
addpath('~/data/code_resources/m_common/linreg/');
%addpath('/home/greg/data/code_resources/m_common/hist2/');

% Set data path and file name, read in file
rawdatapath = '../rawdata/soilsensors_hourly/';
processeddatapath = '../processed_data/';

% LOAD the data (can switch between daily/hourly data here)
climData = csvread([processeddatapath 'wyear_climatesummary.txt'], 1,0);
tsData = csvread([processeddatapath 'wyear_soiltempsummary_hourly.txt'], 1,0);
% tsData = csvread([processeddatapath 'wyear_soiltempsummary_daily.txt']);
% Soil water content data
if strcmpi(normalize, 'y')
    vwcData = csvread([processeddatapath ...
        'wyear_soilwatersummary_hourly_smnorm.txt'], 1,0);
%     vwcData = csvread([processeddatapath ...
%         'wyear_soilwatersummary_daily_smnorm.txt'], 1,0);
else
    vwcData = csvread([processeddatapath...
        'wyear_soilwatersummary_hourly.txt'], 1,0);
%     vwcData = csvread([processeddatapath...
%         'wyear_soilwatersummary_daily.txt'], 1,0);
end

% climData includes more than just soil sites, 
% Get a subset corresponding to the sites and years in tsData
matchsoil = ismember(climData(:, 1:2), tsData(:, 1:2), 'rows');

% Aggregation index for climData 
[allsites, ~, valindex] = unique(climData(:,1));
climAggindex = [valindex ones(size(climData(:,1)))];

% Aggregation index for climData(matchsoil) data
[soilsites, ~, valindex] = unique(climData(matchsoil,1));
soilAggindex = [valindex ones(size(climData(matchsoil,1)))];

% Unique elevations in climData
elevAllAgg = accumarray(climAggindex, climData(:,89), [numel(allsites) 1], @mean);
%Unique elevations in climData(matchsoil)
elevSoilAgg = accumarray(soilAggindex, climData(matchsoil,89), [numel(soilsites) 1], @mean);

% Assign climData variables using the headers file
fid = fopen([processeddatapath 'headersClim.txt']);
headerCell = textscan(fid, '%s', 'headerlines', 1);
fclose(fid);
headers = headerCell{1};
for i=1:11
    eval([headers{i} ' = climData(:,i);']);
end
maxswe = maxswe*25.4;
% Load precip + SWE and convert to mm
for i=12:54
    eval([headers{i} ' = climData(:,i)*25.4;']);
end
% and the rest with no conversion
for i=55:length(headers)
    eval([headers{i} ' = climData(:,i);']);
end

% Assign vwcData variables using the headers file
fid = fopen([processeddatapath 'headersVWC.txt']);
headerCell = textscan(fid, '%s', 'headerlines', 1);
fclose(fid);
headers = headerCell{1};
for i=1:length(headers)
    eval([headers{i} ' = vwcData(:,i);']);
end

% Assign tsData variables using the headers file
fid = fopen([processeddatapath 'headersTsoil.txt']);
headerCell = textscan(fid, '%s', 'headerlines', 1);
fclose(fid);
headers = headerCell{1};

for i=1:length(headers)
    eval([headers{i} ' = tsData(:,i);']);
end


%------------------------------------------------------------------
% FIG 1 - Plot data for entire network and soil subset

% Get the soil site subset of the climate data - note that this includes
% all climate years - so it is different than climData(matchsoil) dataset
testsoil = ismember(siteClim, soilsites);
% Then select years greater than 2000
testyears = yearClim > 2000;

fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Wateryear data scatter - all sites/years');

subplot (4, 2, 1)
plot(elevAllAgg, elevAllAgg, 'ok', 'MarkerFaceColor', 'k');
hold on;
plot(elevSoilAgg, elevSoilAgg+100, 'or');
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
testsoil = ismember(siteClim, soilsites);
plot(elev(testsoil), accumprecip(testsoil), 'or');
ylabel('Annual precip (mm)');
text(0.6, 0.7, 'Annual precip', 'Units', 'normalized');
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
title('Snowpack start day');

subplot (4, 2, 8)
plot(elev, meltdoy, 'ok', 'MarkerFaceColor', 'k');
hold on;
plot(elev(testsoil), meltdoy(testsoil), 'or');
xlabel('Elevation (m)'); ylabel('Day of year');
title('Snow-free day');

%------------------------------------------------------------------
% FIG 2 - Histograms of entire network and soil subset - 2001-2011
% testsoil & testyears is defined in FIG 1
fignum = fignum+1;
fig = figure('position',[100 0 1200 800],'paperpositionmode',...
    'auto', 'color','white','InvertHardcopy','off');
set(fig, 'Name', 'Wateryear data histograms - all sites/years');
set(fig, 'DefaultAxesFontSize',14, 'DefaultTextFontSize', 15);

subplot (4, 2, 1)
xedges = linspace(500, 4000, 60);
networkhist = histc(elevAllAgg, xedges);
soilhist = histc(elevSoilAgg, xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'FaceColor', [0.7 0.7 0.7]);
set(gca, 'position', [0.95 0.975 1.1 1.05] .* get(gca, 'position'));
xlim([700 3700]); ylim([0 40]);
%vline(nanmean(elevAllAgg), '--k');
%h = vline(nanmean(elevSoilAgg), '--'); set(h, 'Color', [0.5 0.5 0.5]);
text(0.05, 0.4, 'Elevation (m)', 'Units', 'normalized');
text(0.95, 0.85, 'a', 'Units', 'normalized');
%ylabel('Frequency');

subplot (4, 2, 2)
xedges = linspace(-5, 20, 60);
networkhist = histc(maat(testyears), xedges);
soilhist = histc(maat(testsoil & testyears), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'FaceColor', [0.7 0.7 0.7]);
set(gca, 'position', [0.95 0.975 1.1 1.05] .* get(gca, 'position'));
%vline(nanmean(maat), '--k');
%h=vline(nanmean(maat(testsoil)), '--'); set(h, 'Color', [0.5 0.5 0.5]);
xlim([-5 15]); ylim([0 500]);
text(0.55, 0.4, 'Mean annual T_{air} (^oC)', 'Units', 'normalized');
text(0.95, 0.85, 'b', 'Units', 'normalized');

subplot (4, 2, 3)
xedges = linspace(0, 2500, 60);
networkhist = histc(accumprecip(testyears), xedges);
soilhist = histc(accumprecip(testsoil & testyears), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'FaceColor', [0.7 0.7 0.7]);
set(gca, 'position', [0.95 0.975 1.1 1.05] .* get(gca, 'position'));
%vline(nanmean(accumprecip), '--k');
%h=vline(nanmean(accumprecip(testsoil)), '--'); set(h, 'Color', [0.5 0.5 0.5]);
xlim([-10 2500]); ylim([0 500]);
text(0.55, 0.4, 'Annual precip (mm)', 'Units', 'normalized');
%ylabel('Frequency');
text(-0.15, -1,'Frequency of ocurrence', 'Units', 'normalized', 'Rotation', 90);
text(0.95, 0.85, 'c', 'Units', 'normalized');

subplot (4, 2, 4)
xedges = linspace(0, 800, 60);
networkhist = histc(JASprecip(testyears), xedges);
soilhist = histc(JASprecip(testsoil & testyears), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'FaceColor', [0.7 0.7 0.7]);
set(gca, 'position', [0.95 0.975 1.1 1.05] .* get(gca, 'position'));
% vline(nanmean(JASprecip), '--k');
% h=vline(nanmean(JASprecip(testsoil)), '--'); set(h, 'Color', [0.5 0.5 0.5]);
xlim([-20 570]); ylim([0 600]);
text(0.535, 0.5, ['Jul+Aug+Sep' 10  'precip (mm)'],...
    'Units', 'normalized');
text(0.95, 0.85, 'd', 'Units', 'normalized');

subplot (4, 2, 5)
xedges = linspace(10, 2000, 60);
networkhist = histc(maxswe(testyears), xedges);
soilhist = histc(maxswe(testsoil & testyears), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'FaceColor', [0.7 0.7 0.7]);
set(gca, 'position', [0.95 0.975 1.1 1.05] .* get(gca, 'position'));
% vline(nanmean(maxswe), '--k');
% h=vline(nanmean(maxswe(testsoil)), '--'); set(h, 'Color', [0.5 0.5 0.5]);
xlim([-2 2000]); ylim([0 400]);
text(0.55, 0.4, 'Peak SWE (mm)', 'Units', 'normalized');
%ylabel('Frequency');
text(0.95, 0.85, 'e', 'Units', 'normalized');

subplot (4, 2, 6)
xedges = linspace(0, 365, 60);
networkhist = histc(totaldaysSC(testyears), xedges);
soilhist = histc(totaldaysSC(testsoil & testyears), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'FaceColor', [0.7 0.7 0.7]);
set(gca, 'position', [0.95 0.975 1.1 1.05] .* get(gca, 'position'));
%vline(nanmean(totaldaysSC), '--k');
%h=vline(nanmean(totaldaysSC(testsoil)), '--'); set(h, 'Color', [0.5 0.5 0.5]);
xlim([-2 300]); ylim([0 500]);
text(0.05, 0.4, 'Total snow-covered days', 'Units', 'normalized');
text(0.95, 0.85, 'f', 'Units', 'normalized');

ticklocs = [0, 32, 62, 93, 121];
tickmonths = ['Oct 1 '; 'Nov 1 '; 'Dec 1 '; 'Jan 1 '; 'Feb 1 '];
subplot (4, 2, 7)
xedges = linspace(0, 130, 60);
networkhist = histc(onsetdoy(testyears), xedges);
soilhist = histc(onsetdoy(testsoil & testyears), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'FaceColor', [0.7 0.7 0.7]);
set(gca, 'position', [0.95 0.975 1.1 1.05] .* get(gca, 'position'));
%vline(nanmean(onsetdoy), '--k');
%h=vline(nanmean(onsetdoy(testsoil)), '--'); set(h, 'Color', [0.5 0.5 0.5]);
xlim([-2 125]); ylim([0 700]);
set(gca,'XTick', ticklocs, 'XTickLabel', tickmonths);
text(0.55, 0.4, 'Snowpack start day', 'Units', 'normalized');
%ylabel('Frequency');
text(0.95, 0.85, 'g', 'Units', 'normalized');

ticklocs = [122, 150, 181, 211, 242, 272, 303, 334];
tickmonths = ['Feb 1 '; 'Mar 1 '; 'Apr 1 '; 'May 1 '; 'Jun 1 ';...
    'Jul 1 '; 'Aug 1 '; 'Sep 1 '];
subplot (4, 2, 8)
xedges = linspace(0, 365, 60);
networkhist = histc(meltdoy(testyears), xedges);
soilhist = histc(meltdoy(testsoil & testyears), xedges);
bar(xedges, networkhist, 'k');
hold on;
bar(xedges, soilhist, 'FaceColor', [0.7 0.7 0.7]);
set(gca, 'position', [0.95 0.975 1.1 1.05] .* get(gca, 'position'));
%vline(nanmean(meltdoy), '--k');
%h=vline(nanmean(meltdoy(testsoil)), '--'); set(h, 'Color', [0.5 0.5 0.5]);
xlim([100 340]); ylim([0 700]);
set(gca,'XTick', ticklocs, 'XTickLabel', tickmonths);
text(0.05, 0.4, 'Snow-free day', 'Units', 'normalized');
text(0.95, 0.85, 'h', 'Units', 'normalized');

figpath = '../figures/';
print(fig,'-depsc2','-painters',[figpath 'figA.eps']) 

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
% testsoil = ismember(siteClim, soilsites);
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
subplot (2, 2, 1);
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
% FIG 4 - Plot MAST vs meltdoy and onsetdoy days by depth and elevation
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'MAST vs MAT/Snowcovered days by depth and elevation');

% Left side - Plot MAST vs onset/melt days by depth
subplot (2, 2, 1)
plot(onsetdoy(matchsoil), mast5cm, '.g', ...
    onsetdoy(matchsoil), mast20cm, '.b',...
    onsetdoy(matchsoil), mast50cm, '.k');
legend('5cm', '20cm', '50cm');
xlabel('Snowpack onset day');
ylabel('MAST(^oC)');
title('Mean wateryear Ts vs Mean wateryear Tair');

subplot (2, 2, 3)
plot(meltdoy(matchsoil), mast5cm, '.g', ...
    meltdoy(matchsoil), mast20cm, '.b',...
    meltdoy(matchsoil), mast50cm, '.k');
xlabel('Snowmelt day');
ylabel('MAST (^oC)');
title('... vs Snowcovered days');

% Right side - Plot MAST vs MAT/Snowcovered days in elevation bins
% First get the elevation categories
testhi = elev(matchsoil) > 3000;
testmid = elev(matchsoil) > 2500 & elev(matchsoil) < 3000;
testlo = elev(matchsoil) > 2000 & elev(matchsoil) < 2500;
% Match the climatedata variables with soil data
matchMelt = meltdoy(matchsoil);
matchOnset = onsetdoy(matchsoil);

subplot (2, 2, 2)
plot(matchOnset(testhi), mast20cm(testhi), '.b', ...
    matchOnset(testmid), mast20cm(testmid), '.g',...
    matchOnset(testlo), mast20cm(testlo), '.r');
legend('3000+', '2500-3000', '2000-2500cm');
xlabel('Snowpack onset day');
ylabel('20cm MAST (^oC)');
title('Mean wateryear Ts vs. Mean wateryear Tair');

subplot (2, 2, 4)
plot(matchMelt(testhi), mast20cm(testhi), '.b', ...
    matchMelt(testmid), mast20cm(testmid), '.g',...
    matchMelt(testlo), mast20cm(testlo), '.r');
xlabel('Snowmelt day');
ylabel('20cm MAST(^oC)');
title('... vs Snowcovered days');

%--------------------------------------------------------------
% FIG 6 - Plot OFFSETS between MAST and MAT
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

subplot (2, 2, 4)
plot(onsetdoy(matchsoil), offset5cm, '.r', ...
    onsetdoy(matchsoil), offset20cm, '.g',...
    onsetdoy(matchsoil), offset50cm, '.b');
hold on;
xrange = [0, 200]; %xlim(gca);
[~, rsq, xfit, yfit] = fitline(onsetdoy(matchsoil), offset20cm, 1, xrange);
plot(xfit, yfit, '--k');
text(15, 7, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('Onset day'); ylabel('^oC');
title('MAST - MAT vs Onset day');
%--------------------------------------------------------------
% FIG 7 - Plot MAST vs snowpack duration for 3 sites
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Regress mast vs snowpack duration for 3 sites');
test = siteClim==828;
site1dur = snowduration(matchsoil & test);
site1maat = maat(matchsoil & test);
site1mast = mast20cm(siteTsoil==828);
test = siteClim==393;
site2dur = snowduration(matchsoil & test);
site2maat = maat(matchsoil & test);
site2mast = mast20cm(siteTsoil==393);
test = siteClim==582;
site3dur = snowduration(matchsoil & test);
site3mast = mast20cm(siteTsoil==582);
site3maat = maat(matchsoil & test);

subplot (3, 2, 1)
plot(site1dur, site1mast, '.b');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site1dur, site1mast, 1, xrange);
plot(xfit, yfit,':k');
text(mean(get(gca,'Xlim')), mean(get(gca,'Ylim')), ...
    ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('Snowpack duration (days)'); ylabel('MAST');
title('Trial Lake');

subplot (3, 2, 3)
plot(site2dur, site2mast, '.b');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site2dur, site2mast, 1, xrange);
plot(xfit, yfit,':k');
text(mean(get(gca,'Xlim')), mean(get(gca,'Ylim')), ...
    ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('Snowpack duration (days)'); ylabel('MAST');
title('Chalk Creek 1');

subplot (3, 2, 5)
plot(site3dur, site3mast, '.b');
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site3dur, site3mast, 1, xrange);
plot(xfit, yfit,':k');
text(mean(get(gca,'Xlim')), mean(get(gca,'Ylim')), ...
    ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('Snowpack duration (days)'); ylabel('MAST');
title('Little Bear');

subplot (3, 2, 2)
plot(site1maat, site1mast, '.b')
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site1maat, site1mast, 1, xrange);
plot(xfit, yfit,':k');
text(mean(get(gca,'Xlim')), mean(get(gca,'Ylim')), ...
    ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('MAT'); ylabel('MAST');
title('Trial Lake');

subplot (3, 2, 4)
plot(site2maat, site2mast, '.b')
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site2maat, site2mast, 1, xrange);
plot(xfit, yfit,':k');
text(mean(get(gca,'Xlim')), mean(get(gca,'Ylim')), ...
    ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('MAT'); ylabel('MAST');
title('Chalk Creek 1');

subplot (3, 2, 6)
plot(site3maat, site3mast, '.b')
hold on;
xrange = xlim(gca);
[~, rsq, xfit, yfit] = fitline(site3maat, site3mast, 1, xrange);
plot(xfit, yfit,':k');
text(mean(get(gca,'Xlim')), mean(get(gca,'Ylim')), ...
    ['r^2 = ' num2str(rsq, 2)]); % r^2 values
xlabel('MAT'); ylabel('MAST');
title('Little Bear');

%----------------------------------------------------
% FIG 8 - Plot snowcovered temp vs onset temps
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
% FIG 9 - Plot soil moisture vs max swe/snowmelt day
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





