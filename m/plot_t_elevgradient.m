% plot_t_elevgradient.m
%
% Plots elevational gradients in air and soil temperature, and the offset
% between the two, including comparisons of growing season vs winter months
% (user selected months used for comparison).
%
% User input selects a subset of sites (by state, or all) and a comparison
% month - will be compared with july.
%


clear;          % clear memory
close all;      % clear any figures
fignum=0;       % used to increment figure number for plots

% Add needed tools
addpath('/home/greg/data/code_resources/m_common/');
addpath('~/data/code_resources/m_common/linreg/');
addpath('/home/greg/data/code_resources/m_common/nanstats/');
addpath('/home/greg/data/code_resources/m_common/hline_vline/');

% Set data path and file name, read in file
rawdatapath = '../rawdata/soilsensors_hourly/';
processeddatapath = '../processed_data/';

% Get an inventory of siteIDs and their states from inventory file
% Create format string (station,elev,cdbs_id only here)
formatstr = '%*s%*s%*s%*s%s%f%f%f%f';
fid = fopen([processeddatapath 'SNOTELinventory.csv']);
inventorycell = textscan(fid, formatstr,'Headerlines', 1, 'Delimiter', ',');
fclose(fid);

states = inventorycell{1};
siteIDs = inventorycell{2};
clear inventorycell;

% LOAD the data (can switch between daily/hourly data here)
climData = csvread([processeddatapath 'wyear_climatesummary.txt'], 1,0);
tsData = csvread([processeddatapath 'wyear_soiltempsummary_hourly.txt'], 1,0);
% tsData = csvread([processeddatapath 'wyear_soiltempsummary_daily.txt']);

% climData includes more than just soil sites, 
% Get a subset corresponding to the sites and years in tsData
matchsoil = ismember(climData(:, 1:2), tsData(:, 1:2), 'rows');

% Aggregation index for climData 
[soilsites, ~, valindex] = unique(climData(matchsoil,1));
aggindex = [valindex ones(size(valindex))];
%aggindex = [valindex ones(size(climData(matchsoil,1)))];

% Assign climData variables using the headers file - USE MATCHSOIL
fid = fopen([processeddatapath 'headersClim.txt']);
headerCell = textscan(fid, '%s', 'headerlines', 1);
fclose(fid);
headers = headerCell{1};
for i=1:11
    eval([headers{i} ' = climData(matchsoil,i);']);
end
maxswe = maxswe*25.4;
% Load precip + SWE and convert to mm
for i=12:54
    eval([headers{i} ' = climData(matchsoil,i)*25.4;']);
end
% and the rest with no conversion
for i=55:length(headers)
    eval([headers{i} ' = climData(matchsoil,i);']);
end

% Assign tsData variables using the headers file
fid = fopen([processeddatapath 'headersTsoil.txt']);
headerCell = textscan(fid, '%s', 'headerlines', 1);
fclose(fid);
headers = headerCell{1};
for i=1:length(headers)
    eval([headers{i} ' = tsData(:,i);']);
end

% Aggregate a bunch of variables for plotting
siteClimAgg = accumarray(aggindex, siteClim, [numel(soilsites) 1], @mean);
elevAgg = accumarray(aggindex, elev, [numel(soilsites) 1], @mean);
% Tair and StdDev
julTairMeanAgg = accumarray(aggindex, julTairMean, [numel(soilsites) 1], @nanmean);
julTairSdAgg = accumarray(aggindex, julTairSd, [numel(soilsites) 1], @mean);
janTairMeanAgg = accumarray(aggindex, janTairMean, [numel(soilsites) 1], @nanmean);
janTairSdAgg = accumarray(aggindex, janTairSd, [numel(soilsites) 1], @mean);
%Tsoil (mean and StdDev) for 5 and 20cm
julTs5meanAgg = accumarray(aggindex, julTs5mean, [numel(soilsites) 1], @nanmean);
julTs20meanAgg = accumarray(aggindex, julTs20mean, [numel(soilsites) 1], @nanmean);
julTs5sdAgg = accumarray(aggindex, julTs5sd, [numel(soilsites) 1], @mean);
julTs20sdAgg = accumarray(aggindex, julTs20sd, [numel(soilsites) 1], @mean);
janTs5meanAgg = accumarray(aggindex, janTs5mean, [numel(soilsites) 1], @nanmean);
janTs20meanAgg = accumarray(aggindex, janTs20mean, [numel(soilsites) 1], @nanmean);
janTs5sdAgg = accumarray(aggindex, janTs5sd, [numel(soilsites) 1], @mean);
janTs20sdAgg = accumarray(aggindex, janTs20sd, [numel(soilsites) 1], @mean);

% PLOTS
%----------------------------------
% FIG 1 - January and July Tair/Tsoil gradients - All SNOTEL/years
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name','Jan/July Air and soil temperature gradients, All SNOTEL');

subplot (2, 2, 1)
plot(elev, julTairMean, 'ok', 'MarkerFaceColor', 'r');
hold on
plot(elev, julTs20mean, 'ok', 'MarkerFaceColor', 'k');
xlabel('Elevation(m)');
ylabel('Mean July temp (^oC)');
title('July');
legend('Air', 'Soil(20cm)');

subplot (2, 2, 2)
plot(elev, janTairMean, 'ok', 'MarkerFaceColor', 'b');
hold on
plot(elev, janTs20mean, 'ok', 'MarkerFaceColor', 'w');
% Plot moist adiabatic lapse rate
plot([900 3500], [4, -9], ':k')
xlabel('Elevation(m)');
ylabel('Mean temp (Celsius)');
title('January');
legend('Air', 'Soil', 'Moist adiabatic lapse( 5^oC/km)');

subplot (2, 2, 3)
plot(elevAgg, julTairMeanAgg, 'ok', 'MarkerFaceColor', 'r');
hold on
plot(elevAgg, julTs20meanAgg, 'ok', 'MarkerFaceColor', 'k');
xlabel('Elevation(m)');
ylabel('Mean July temp (^oC)');
title('July (mean of all wateryears)');
legend('Air', 'Soil(20cm)');

subplot (2, 2, 4)
plot(elevAgg, janTairMeanAgg, 'ok', 'MarkerFaceColor', 'b');
hold on
plot(elevAgg, janTs20meanAgg, 'ok', 'MarkerFaceColor', 'w');
% Plot moist adiabatic lapse rate
plot([900 3500], [4, -9], ':k')
xlabel('Elevation(m)');
ylabel('Mean temp (Celsius)');
title('January (mean of all wateryears)');
legend('Air', 'Soil', 'Moist adiabatic lapse( 5^oC/km)');

% -------------------------------------------------------------
% FIG 2
% % Plot all temperature data
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Soil temperature elevation/swe gradients - All SNOTEL/years');

% Plot Snowcovered and snow-free soil temps
subplot(2,2,1);
errorbar(elev, snowcovTs20mean, snowcovTs20sd, 'b.');
hold on
errorbar(elev, snowfreeTs20mean, snowfreeTs20sd, 'r.');
title('Snowcovered and snowfree temperature');
legend('Snow-covered', 'Snow-free');
ylabel('20cm soil temp');
xlabel('Elevation (m)');

subplot(2,2,2);
errorbar(maxswe, snowcovTs20mean, snowcovTs20sd, 'b.');
hold on
errorbar(maxswe, snowfreeTs20mean, snowfreeTs20sd, 'r.');
title('Snowcovered and snowfree temperature');
legend('Snow-covered', 'Snow-free');
ylabel('20cm Ts (^oC)');
xlabel('Peak SWE (mm)');

% Plot January and July soil temps by elevation
subplot(2,2,3);
errorbar(elev, janTs20mean,janTs20sd, 'ko', 'MarkerFaceColor', 'w');
hold on;
errorbar(elev, julTs20mean, julTs20sd, 'ko', 'MarkerFaceColor', 'k');
% Plot moist adiabatic lapse rate
plot([1600 3400], [18, 9], ':k');
legend('January', 'July', 'Moist adiabatic lapse (5^oC/km)');
xlabel('Elevation (m)');
ylabel('20cm Ts (^oC)')

% Plot January and July soil temps by elevation
subplot(2,2,4);
errorbar(maxswe, janTs20mean,janTs20sd, 'ko', 'MarkerFaceColor', 'w');
hold on;
errorbar(maxswe, julTs20mean, julTs20sd, 'ko', 'MarkerFaceColor', 'k');
% Plot moist adiabatic lapse rate
legend('January', 'July', 'Moist adiabatic lapse (5^oC/km)');
xlabel('Peak SWE (mm)');
ylabel('20cm Ts (^oC)')

%----------------------------------
% FIG 3 - Other landscape gradients
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name', 'Other landscape gradients');

subplot (2,2,1)
plot(elev, maat, '.', 'Color', [0.7 0.7 0.7]);
hold on
plot(elevAgg, julTairMeanAgg, '.m');
plot(elevAgg, janTairMeanAgg, '.b');
xlabel('Elevation(m)');
ylabel('Mean T (^oC)');
xlim([700 3700]);
title('Temperature-elevation gradients');
legend('Wateryear MAT', 'Jul MAT', 'Jan MAT');

subplot (2,2,2)
plot(elev, maxswe, '.', 'Color', ...
    [0.7 0.7 0.7]);
hold on
plot(elev, ltMeanSWE, 'ok', 'MarkerFaceColor', 'w');
plot(elev, ltMeanPrecip, 'ok', 'MarkerFaceColor', 'b');
xlabel('Elevation(m)');
ylabel('mm H_2O');
xlim([700 3700]);
title('Precip/SWE-elevation gradients');
legend('PeakSWE', 'LT Peak SWE', 'LT Precip');

subplot(2,2,3)
% Aggregate latitude
latAgg = accumarray(aggindex, lat, [numel(soilsites) 1], @mean);
plot(lat, maat, '.', 'Color', ...
    [0.7 0.7 0.7]);
hold on
plot(latAgg, julTairMeanAgg, '.', 'Color', [0.4 0.4 0.4]);
plot(latAgg, janTairMeanAgg, '.k');
xlabel('Latitude');
ylabel('Mean T (^oC)');
title('Mean temperature vs latitude');
legend('Wateryear MAT', 'Jul MAT', 'Jan MAT');

subplot (2,2,4)
plot(lat, onsetdoy, '.', 'Color', ...
    [0.7 0.7 0.7]);
hold on
plot(lat, meltdoy, '.', 'Color', [0.4 0.4 0.4]);
plot(lat, totaldaysSC, '.k');
xlabel('Latitude');
ylabel('Day of wateryear');
title('Snow events vs latitude');
legend('Snow onset', 'Snow melt', 'Total snow days');

% -----------------------------------------------------------
% FIG 4 - Month soil temps vs air temp
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name', 'Mean monthly air and soil T - All SNOTEL/years');

% Set some plotting parameters
plotorder = 1:3;
months = ['oct';'nov';'dec';'jan';'feb';'mar';'apr';'may';...
    'jun';'jul';'aug';'sep'] ;
polyorder = 1;

for i=plotorder
    subplot(4, 3, i);
    x = elev;
    eval(['y1 = ' months(i,:) 'TairMean;']);
    eval(['y2 = ' months(i,:) 'Ts5mean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
%     [~, rsq, xfit, yfit] = fitline(x, y1, polyorder, [0, 1500]);
%     plot(xfit, yfit, '-k', 'LineWidth', 2);
%     text(700, 8, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
    plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    title(months(i,:));
    xlim([800, 4000]); ylim([-15, 30]);
    hline(0, '--k');
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('^oC')
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
    
    subplot(4,3, i+3)
    eval(['y1 = ' months(i+3,:) 'TairMean;']);
    eval(['y2 = ' months(i+3,:) 'Ts20mean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    title(months(i+3,:));
    xlim([800, 4000]); ylim([-15, 30]);
    hline(0, '--k');
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('^oC')
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
    
    subplot(4,3, i+6);
    eval(['y1 = ' months(i+6,:) 'TairMean;']);
    eval(['y2 = ' months(i+6,:) 'Ts50mean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    title(months(i+6,:));
    xlim([800, 4000]); ylim([-15, 30]);
    hline(0, '--k');
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('^oC')
    elseif i>1
        set(gca, 'YTickLabel', '');
    end

    subplot(4,3, i+9);
    eval(['y1 = ' months(i+9,:) 'TairMean;']);
    eval(['y2 = ' months(i+9,:) 'Ts50mean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    title(months(i+9,:));
    xlim([800, 4000]); ylim([-15, 30]);
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

% -------------------------------------------------------
% FIG 5 - Month soil temps vs air temp
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name', 'Dec-Apr air and soil temp - All SNOTEL/years');

% Set some plotting parameters
plotorder = 1:5;
months = ['dec';'jan';'feb';'mar';'apr'] ;
polyorder = 1;

for i=plotorder
    subplot(3, 5, i);
    x = elev;
    eval(['y1 = ' months(i,:) 'TairMean;']);
    eval(['y2 = ' months(i,:) 'Ts5mean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
%     [~, rsq, xfit, yfit] = fitline(x, y1, polyorder, [0, 1500]);
%     plot(xfit, yfit, '-k', 'LineWidth', 2);
%     text(700, 8, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
    plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    xlim([800, 4000]); ylim([-10, 10]);
    hline(0, '--k');
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('^oC')
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
    
    subplot(3, 5, i+5)
    eval(['y2 = ' months(i,:) 'Ts20mean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    xlim([800, 4000]); ylim([-10, 10]);
    hline(0, '--k');
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('^oC')
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
    
    subplot(3, 5, i+10);
    eval(['y2 = ' months(i,:) 'Ts50mean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    xlim([800, 4000]); ylim([-10, 10]);
    hline(0, '--k');
    if i==1
        ylabel('^oC')
    elseif i==3
        xlabel('Elevation (m)');
        set(gca, 'YTickLabel', '');
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
end

%--------------------------------------------------------
% FIG 6 - Restricted gradients
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name','Restricted gradients');

%Select peak SWE gradient
swetest = ltMeanSWE<600 & ltMeanSWE>400;

subplot (2, 2, 1)
plot(ltMeanSWE(swetest), julTairMean(swetest), '.k');
xlabel('Mean peak SWE(mm)');
ylabel('Mean July temp (Celsius)');
title('Constant SWE - July air temp');
%legend('Air', 'Soil(20cm, 12am)');

subplot (2, 2, 2)
plot(ltMeanSWE(swetest), elev(swetest), '.k');
xlabel('Mean peak SWE(mm)');
ylabel('elevation(m)');
title('Constant SWE - elevation');

% Select an elevation/temp gradient
temptest = elev<2550 & elev>2350;

subplot (2, 2, 3)
plot(elev(temptest), julTairMean(temptest), '.k');
xlabel('elevation(m)');
ylabel('Mean July temperature (Celsius');
title('Constant elevation - July temperature');

subplot (2, 2, 4)
plot(elev(temptest), ltMeanSWE(temptest), '.k');
xlabel('elevation(m)');
ylabel('Mean peak SWE(mm)');
title('Constant elevation - Peak SWE');

%----------------------------------------------------------------------
% FIG 7 - Winter Air/Soil gradient - All SNOTEL/years
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name','Dec Air/Soil gradient - All SNOTEL/years');

subplot (1,2,1);
plot(elev, decTairMean, 'ok', 'MarkerFaceColor', [0.7 0.7 0.7]);
hold on
plot(elev, decTs20mean, 'ok', 'MarkerFaceColor', 'w');
% Plot moist adiabatic lapse rate
plot([900 3500], [4, -9], ':k')
xlabel('Elevation(m)');
ylabel('Mean temp (Celsius)');
legend('Air', 'Soil(20cm)', 'Moist adiabatic lapse( 5^oC/km)');

subplot (1,2,2);
plot(elev, decTairMean, 'ok', 'MarkerFaceColor', [0.7 0.7 0.7]);
hold on
plot(elev, decTs5mean, 'ok', 'MarkerFaceColor', 'w');
% Plot moist adiabatic lapse rate
plot([900 3500], [4, -9], ':k')
xlabel('Elevation(m)');
ylabel('Mean temp (Celsius)');
legend('Air', 'Soil(5cm)', 'Moist adiabatic lapse( 5^oC/km)');

% -------------------------------------------------------------
% FIG 8 - Winter Air/Soil gradient and offset - All SNOTEL/years
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name','Dec Air/Soil gradient and offset - All SNOTEL/years');

subplot (1,2,1);
errorbar(elev, decTairMean, decTairSd, 'ok',...
    'MarkerFaceColor', [0.7 0.7 0.7]);
hold on
errorbar(elev, decTs20mean, decTs20sd, 'ok', 'MarkerFaceColor', 'w');
% Plot moist adiabatic lapse rate
plot([900 3500], [4, -9], ':k')
xlabel('Elevation(m)');
ylabel('Mean temp (^oC)');
legend('Air', 'Soil (20cm)', 'Moist adiabatic lapse( 5^oC/km)');
title('Air & soil temp.');

subplot(1,2,2);
offset = abs(decTairMean - decTs20mean);
plot(elev, offset(:,1), 'ok', 'MarkerFaceColor', 'b');
hold on;

xfit = linspace(900, 3500);
[b,bint,resid,rint,stats] = regress2(offset, [elev ones(size(elev))]);
handles(2) = plot(xfit, polyval(b, xfit), '--k', 'Linewidth', 1.5);
%text(1750, -7, ['y = ' num2str(b(1),'%1.4f') 'x + ' num2str(b(2),'%2.1f')]);
if stats(3) < 0.01
    text(1750, 2, ['r^2 = ' num2str(stats(1),2) ', p < 0.01']);
else
    text(1750, 2, ['r^2 = ' num2str(stats(1),2)...
        ', p = ' num2str(stats(3),2)]);
end

xlabel('Elevation (m)');
ylabel('^oC');
title('Air-Soil temp. offset');
%legend('Air-Soil Offset');

% -------------------------------------------------------------
% FIG 9 - Jan Air/Soil T gradients - Utah data only
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name','(UT) January Air/Soil T gradients');

% Get selected state site out of all data
selectState = 'CO';
selectIDs = siteIDs(strcmpi(states, selectState));
st_test = ismember(siteClim, selectIDs);

subplot (1,1,1);
errorbar(elev(st_test), janTairMean(st_test), janTairSd(st_test), ...
    'ok', 'MarkerFaceColor', [0.7 0.7 0.7], 'Color', [0.7 0.7 0.7], ...
    'MarkerEdgeColor', 'k');
hold on
errorbar(elev(st_test), janTs20mean(st_test), janTs20sd(st_test), ...
    'ok', 'MarkerFaceColor', 'w', 'Color', [0.7 0.7 0.7], ...
    'MarkerEdgeColor', 'k');
% Plot moist adiabatic lapse rate
plot([900 3500], [4, -9], ':k')
xlim([1500 3500]);
xlabel('Elevation(m)');
ylabel('Mean temp (^oC)');
legend('Air', 'Soil (20cm)', 'Moist adiabatic lapse( 5^oC/km)');
title('Utah - January Tair/Tsoil');

% -------------------------------------------------------------
% FIG 10 - Jan Air/Soil T gradients - aggregated Utah data only
fignum = fignum+1;    
h = figure('position',[100 0 650 600],'paperpositionmode',...
    'auto', 'color','white','InvertHardcopy','off');
set(h, 'Name','(UT-Agg) January Air/Soil T gradients',...
    'DefaultAxesFontSize',18, 'DefaultTextFontSize', 18);

% Get selected state site out of aggregated data
st_testAgg = ismember(siteClimAgg, selectIDs);

% Assign to x and y variables
x = elevAgg(st_testAgg);
y1 = janTairMeanAgg(st_testAgg);
sd1 = janTairSdAgg(st_testAgg);
y2 = janTs20meanAgg(st_testAgg);
sd2 = janTs20sdAgg(st_testAgg);

% Mean StDev of Tair
meanStdDev = nanmean(sd1);

% Plot Jan mean Tair and regression line
handles(1) = plot(x, y1, 'ok', 'MarkerFaceColor', [0.7 0.7 0.7],...
    'MarkerEdgeColor', 'k','MarkerSize', 10);
hold on;
xfit = linspace(1692, 3400);
[b,bint,resid,rint,stats] = regress2(y1, [x ones(size(x))]);
handles(2) = plot(xfit, polyval(b, xfit), '--k', 'Linewidth', 1.5);
%text(1750, -7, ['y = ' num2str(b(1),'%1.4f') 'x + ' num2str(b(2),'%2.1f')]);
if stats(3) < 0.01
    text(1750, -7.4, ['r^2 = ' num2str(stats(1),2) ', p < 0.01']);
else
    text(1750, -7.4, ['r^2 = ' num2str(stats(1),2)...
        ', p = ' num2str(stats(3),2)]);
end

% Plot Tair example error bar in an inset box
handles(3) = errorbar(1625, -4.7, meanStdDev,'ok',...
    'MarkerFaceColor',[0.7 0.7 0.7],...
    'MarkerEdgeColor', 'k', 'MarkerSize', 10);
% Create rectangle
annotation(h, 'rectangle', [0.154 0.146 0.05 0.55], 'FaceColor','flat');

% Plot Jan mean Tsoil and regression
handles(4) = errorbar(x, y2, sd2, 'ok', 'MarkerFaceColor', 'White',...
    'MarkerEdgeColor', 'k','MarkerSize', 10);
xfit = linspace(1600, 3400);
[b,bint,resid,rint,stats] = regress2(y2, [x ones(size(x))]);
handles(5) = plot(xfit, polyval(b, xfit), '--k', 'Linewidth', 1.5);
%text(1570, 3.1, ['y = ' num2str(b(1),'%1.4f') 'x + ' num2str(b(2),'%2.1f')]);
if stats(3) < 0.01
    text(1570, 2.2, ['r^2 = ' num2str(stats(1),2) ', p < 0.01']);
else
    text(1570, 2.2, ['r^2 = ' num2str(stats(1),2)...
        ', p = ' num2str(stats(3),2)]);
end

xlim([1500 3500]);
set(gca, 'Ytick', [-8;-6;-4;-2;0;2;4]);
xlabel('Elevation(m)');
ylabel('Mean January T (^oC)');
legend(handles([1 4]), {'T_{air}',...
    'T_{soil}'});
%title('Utah - January Tair/Tsoil (aggregated)');
% Plot moist adiabatic lapse rate
% plot([900 3500], [4, -9], ':k')

figpath = '../figures/';
print(h,'-depsc2','-painters',[figpath 'figD_old.eps']) 

% -------------------------------------------------------------
% FIG 11 - January/July soil T gradients - Utah data only

fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name','(UT) January/July Soil T gradients');
% Plot January and July soil temps by elevation
errorbar(elev(st_test), janTs20mean(st_test),...
    janTs20sd(st_test), 'ko', 'MarkerFaceColor', 'w');
hold on;
errorbar(elev(st_test), julTs20mean(st_test),...
    julTs20sd(st_test), 'ko', 'MarkerFaceColor', 'k');
% Plot moist adiabatic lapse rate
plot([1600 3400], [18, 9], ':k');
legend('January', 'July', 'Moist adiabatic lapse (5^oC/km)');
xlabel('Elevation (m)');
ylabel('20cm Ts (^oC)')

% -------------------------------------------------------------
% FIG 12 - January/July soil T gradients - AggregatedUtah data only

% Assign x and y variables
x = elevAgg(st_testAgg);
y1 = janTs20meanAgg(st_testAgg);
sd1 = janTs20sdAgg(st_testAgg);
y2 = julTs20meanAgg(st_testAgg);
sd2 = julTs20sdAgg(st_testAgg);

fignum = fignum+1;    
h = figure('position',[100 0 650 600],'paperpositionmode',...
    'auto', 'color','white','InvertHardcopy','off');
set(h, 'Name','(UT-Agg) January/July Soil T gradients',...
    'DefaultAxesFontSize',18, 'DefaultTextFontSize', 18);
% Plot January Ts by elevation and regression
handles(1) = errorbar(x, y1, sd1, 'ok', 'MarkerFaceColor', 'White',...
    'MarkerSize', 10);
hold on;
xfit = linspace(1600, 3400);
[b,bint,resid,rint,stats] = regress2(y1, [x ones(size(x))]);
handles(2) = plot(xfit, polyval(b, xfit), '--k', 'Linewidth', 1.5);
%text(1600, 4.7, ['y = ' num2str(b(1),'%1.4f') 'x + ' num2str(b(2),'%2.1f')]);
if stats(3) < 0.01
    text(1600, 3.1, ['r^2 = ' num2str(stats(1),2) ', p < 0.01']);
else
    text(1600, 3.1, ['r^2 = ' num2str(stats(1),2)...
        ', p = ' num2str(stats(3),2)]);
end
% Plot July Ts by elevation and regression
handles(3) = errorbar(x, y2, sd2, 'ok', 'MarkerFaceColor', 'k',...
    'MarkerSize', 10);
[b,bint,resid,rint,stats] = regress2(y2, [x ones(size(x))]);
handles(4) = plot(xfit, polyval(b, xfit), '--k', 'Linewidth', 1.5);
%text(1600, 9.5, ['y = ' num2str(b(1),'%1.4f') 'x + ' num2str(b(2),'%2.1f')]);
if stats(3) < 0.01
    text(1600, 9, ['r^2 = ' num2str(stats(1),2) ', p < 0.01']);
else
    text(1600, 9, ['r^2 = ' num2str(stats(1),2)...
        ', p = ' num2str(stats(3),2)]);
end
% Plot moist adiabatic lapse rate
%plot([1600 3400], [18, 9], ':k');
xlim([1500 3500]);
set(gca, 'Ytick', [0;5;10;15;20;25]);
legend(handles([1 3]), {'January', 'July'});
xlabel('Elevation (m)');
ylabel('Mean monthly T_{soil} (^oC)')

figpath = '../figures/';
print(h,'-depsc2','-painters',[figpath 'figC_old.eps']) 


% -------------------------------------------------------------
% FIG 13 - January/July Tsoil and Tair gradients - AggregatedUtah data only

% Assign x and y variables
x = elevAgg(st_testAgg);
y1 = janTs5meanAgg(st_testAgg);
sd1 = janTs5sdAgg(st_testAgg);
y2 = julTs5meanAgg(st_testAgg);
sd2 = julTs5sdAgg(st_testAgg);

fignum = fignum+1;    
h = figure('position',[100 0 1100 500],'paperpositionmode',...
    'auto', 'color','white','InvertHardcopy','off');
set(h, 'Name','(UT-Agg) January/July Air & Soil T gradients',...
    'DefaultAxesFontSize',18, 'DefaultTextFontSize', 18);

subplot(1,2,1);
% Plot January Ts by elevation and regression
handles(1) = errorbar(x, y1, sd1, 'ok', 'MarkerFaceColor', 'White',...
    'MarkerSize', 10);
hold on;
xfit = linspace(1600, 3400);
[b,bint,resid,rint,stats] = regress2(y1, [x ones(size(x))]);
handles(2) = plot(xfit, polyval(b, xfit), '--k', 'Linewidth', 1.5);
if stats(3) < 0.01
    text(1600, -3.5, ['r^2 = ' num2str(stats(1),2) ', p < 0.01']);
else
    text(1600, -3.5, ['r^2 = ' num2str(stats(1),2)...
        ', p = ' num2str(stats(3),2)]);
end
% Plot July Ts by elevation and regression
handles(3) = errorbar(x, y2, sd2, 'ok', 'MarkerFaceColor', 'k',...
    'MarkerSize', 10);
[b,bint,resid,rint,stats] = regress2(y2, [x ones(size(x))]);
handles(4) = plot(xfit, polyval(b, xfit), '--k', 'Linewidth', 1.5);
if stats(3) < 0.01
    text(1600, 9, ['r^2 = ' num2str(stats(1),2) ', p < 0.01']);
else
    text(1600, 9, ['r^2 = ' num2str(stats(1),2)...
        ', p = ' num2str(stats(3),2)]);
end
set(gca, 'position', [0.90 1 1.15 1] .* get(gca, 'position'));
xlim([1500 3500]);ylim([-13 30]);
set(gca,'Ytick',[-10;-5;0;5;10;15;20;25;30],'Xtick',[1500;2000;2500;3000]);
legend(handles([1 3]), {'January', 'July'}, 'Location', 'Southeast');
text(0.8, 0.9, 'T_{soil}', 'Units', 'normalized', 'Fontangle', 'italic',...
    'Fontsize', 20);
text(0.85, -0.1, 'Elevation (m)', 'Units', 'normalized');
%xlabel('Elevation (m)');
ylabel('Mean monthly T (^oC)')

% Assign x and y variables
y1 = janTairMeanAgg(st_testAgg);
sd1 = janTairSdAgg(st_testAgg);
y2 = julTairMeanAgg(st_testAgg);
sd2 = julTairSdAgg(st_testAgg);

subplot(1,2,2);
% Plot January Tair by elevation and regression
handles(1) = errorbar(x, y1, sd1, 'ok', 'MarkerFaceColor', 'White',...
    'MarkerSize', 10);
hold on;
xfit = linspace(1600, 3400);
[b,bint,resid,rint,stats] = regress2(y1, [x ones(size(x))]);
handles(2) = plot(xfit, polyval(b, xfit), '--k', 'Linewidth', 1.5);
if stats(3) < 0.01
    text(1600, 5.2, ['r^2 = ' num2str(stats(1),2) ', p < 0.01']);
else
    text(1600, 5.2, ['r^2 = ' num2str(stats(1),2)...
        ', p = ' num2str(stats(3),2)]);
end
% Plot July Tair by elevation and regression
handles(3) = errorbar(x, y2, sd2, 'ok', 'MarkerFaceColor', 'k',...
    'MarkerSize', 10);
[b,bint,resid,rint,stats] = regress2(y2, [x ones(size(x))]);
handles(4) = plot(xfit, polyval(b, xfit), '--k', 'Linewidth', 1.5);
if stats(3) < 0.01
    text(1600, 11, ['r^2 = ' num2str(stats(1),2) ', p < 0.01']);
else
    text(1600, 11, ['r^2 = ' num2str(stats(1),2)...
        ', p = ' num2str(stats(3),2)]);
end
% Plot moist adiabatic lapse rate
%plot([1600 3400], [18, 9], ':k');
set(gca, 'position', [0.90 1 1.15 1] .* get(gca, 'position'));
xlim([1500 3500]);ylim([-13 27]);
set(gca, 'YtickLabel', [],'Xtick',[2000;2500;3000;3500]);
text(0.8, 0.9, 'T_{air}', 'Units', 'normalized', 'Fontangle', 'italic',...
    'Fontsize',20);

figpath = '../figures/';
print(h,'-depsc2','-painters',[figpath 'figC_5cm.eps']) 


junk = 99;
