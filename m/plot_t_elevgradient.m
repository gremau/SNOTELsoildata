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
addpath('/home/greg/data/code_resources/m_common/linear/');
addpath('/home/greg/data/code_resources/m_common/nanstuff/');
addpath('/home/greg/data/code_resources/m_common/hline_vline/');

% Set data path and file name, read in file
rawdatapath = '../rawdata/';
processeddatapath = '../processed_data/';

% Load lists of sites with data in the daily/hourly data directory
dailysites = sortrows(csvread('../rawdata/allsensors_daily/filelist.txt'));
soilsites = sortrows(csvread('../rawdata/soilsensors_hourly/filelist.txt'));
sites = unique(soilsites(:,1));

% LOAD the data (can switch between daily/hourly data here)
climData = dlmread([processeddatapath 'wyear_climatesummary.txt'], ',', 1,0);
tsData = csvread([processeddatapath 'wyear_soiltempsummary_hourly.txt'], 1,0);
% tsData = csvread([processeddatapath 'wyear_soiltempsummary_daily.txt']);

% Get a subset of climData that corresponds with available soildata
[matchsoil, idx] = ismember(climData(:, 1:2), tsData(:, 1:2), 'rows');
soilClim = climData(matchsoil, :);
% matchsoil2 = ismember(tsData(:, 1:2), soilClim(:, 1:2), 'rows');

% Ask user for month number and state
% statesel = input(...
%     'Which state ("AZ, CO, ID, MT, NM, NV, UT, WY, or all")?: ', 's');

% Get an inventory of siteIDs and their states from inventory file
% Create format string (station,elev,cdbs_id only here)
formatstr = '%*s%*s%*s%*s%s%f%f%f%f';
fid = fopen([processeddatapath 'SNOTELinventory.csv']);
inventorycell = textscan(fid, formatstr,'Headerlines', 1, 'Delimiter', ',');
fclose(fid);

states = inventorycell{1};
siteIDs = inventorycell{2};
clear inventorycell;

% Aggregation index for soilClim data
[vals, ~, valindex] = unique(soilClim(:,1));
aggindex = [valindex ones(size(soilClim(:,1)))];

%--------------------------------------------------------------
% Assign variables
sites_cl = soilClim(:,1);
maxswe = soilClim(:,3)*25.4;
maat = soilClim(:, 74); % Mean wateryear air temp
onsetdoy = soilClim(:, 6);
meltdoy = soilClim(:, 7);
totaldaysSC = soilClim(:, 9);
octTairMean = soilClim(:,50);
octTairSd = soilClim(:, 51);
novTairMean = soilClim(:, 52);
novTairSd = soilClim(:, 53);
decTairMean = soilClim(:, 54);
decTairSd = soilClim(:, 55);
janTairMean = soilClim(:, 56);
janTairSd = soilClim(:, 57);
febTairMean = soilClim(:, 58);
febTairSd = soilClim(:, 59);
marTairMean = soilClim(:, 60);
marTairSd = soilClim(:, 61);
aprTairMean = soilClim(:, 62);
aprTairSd = soilClim(:, 63);
mayTairMean = soilClim(:, 64);
mayTairSd = soilClim(:, 65);
junTairMean = soilClim(:, 66);
junTairSd = soilClim(:, 67);
julTairMean = soilClim(:, 68);
julTairSd = soilClim(:, 69);
augTairMean = soilClim(:, 70);
augTairSd = soilClim(:, 71);
sepTairMean = soilClim(:, 72);
sepTairSd = soilClim(:, 73);

elev = soilClim(:, 82);
lat = soilClim(:, 83);
lon = soilClim(:, 84);
ltMeanSWE = soilClim(:, 85);
ltMeanPrecip = soilClim(:, 86);

octTs5mean = tsData(:, 3);
octTs5sd = tsData(:, 4);
octTs20mean = tsData(:, 5);
octTs20sd = tsData(:, 6);
octTs50mean = tsData(:, 7);
octTs50sd = tsData(:, 8);
novTs5mean = tsData(:, 9);
novTs5sd = tsData(:, 10);
novTs20mean = tsData(:, 11);
novTs20sd = tsData(:, 12);
novTs50mean = tsData(:, 13);
novTs50sd = tsData(:, 14);
decTs5mean = tsData(:, 15);
decTs5sd = tsData(:, 16);
decTs20mean = tsData(:, 17);
decTs20sd = tsData(:, 18);
decTs50mean = tsData(:, 19);
decTs50sd = tsData(:, 20);
janTs5mean = tsData(:, 21);
janTs5sd = tsData(:, 22);
janTs20mean = tsData(:, 23);
janTs20sd = tsData(:, 24);
janTs50mean = tsData(:, 25);
janTs50sd = tsData(:, 26);
febTs5mean = tsData(:, 27);
febTs5sd = tsData(:, 28);
febTs20mean = tsData(:, 29);
febTs20sd = tsData(:, 30);
febTs50mean = tsData(:, 31);
febTs50sd = tsData(:, 32);
marTs5mean = tsData(:, 33);
marTs5sd = tsData(:, 34);
marTs20mean = tsData(:, 35);
marTs20sd = tsData(:, 36);
marTs50mean = tsData(:, 37);
marTs50sd = tsData(:, 38);
aprTs5mean = tsData(:, 39);
aprTs5sd = tsData(:, 40);
aprTs20mean = tsData(:, 41);
aprTs20sd = tsData(:, 42);
aprTs50mean = tsData(:, 43);
aprTs50sd = tsData(:, 44);
mayTs5mean = tsData(:, 45);
mayTs5sd = tsData(:, 46);
mayTs20mean = tsData(:, 47);
mayTs20sd = tsData(:, 48);
mayTs50mean = tsData(:, 49);
mayTs50sd = tsData(:, 50);
junTs5mean = tsData(:, 51);
junTs5sd = tsData(:, 52);
junTs20mean = tsData(:, 53);
junTs20sd = tsData(:, 54);
junTs50mean = tsData(:, 55);
junTs50sd = tsData(:, 56);
julTs5mean = tsData(:, 57);
julTs5sd = tsData(:, 58);
julTs20mean = tsData(:, 59);
julTs20sd = tsData(:, 60);
julTs50mean = tsData(:, 61);
julTs50sd = tsData(:, 62);
augTs5mean = tsData(:, 63);
augTs5sd = tsData(:, 64);
augTs20mean = tsData(:, 65);
augTs20sd = tsData(:, 66);
augTs50mean = tsData(:, 67);
augTs50sd = tsData(:, 68);
sepTs5mean = tsData(:, 69);
sepTs5sd = tsData(:, 70);
sepTs20mean = tsData(:, 71);
sepTs20sd = tsData(:, 72);
sepTs50mean = tsData(:, 73);
sepTs50sd = tsData(:, 74);

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

% Aggregate by wateryear using accumarray
elevAgg = accumarray(aggindex, elev, [numel(sites) 1], @mean);
julTairMeanAgg = accumarray(aggindex, julTairMean, [numel(sites) 1], @nanmean);
julTs20meanAgg = accumarray(aggindex, julTs20mean, [numel(sites) 1], @nanmean);
janTairMeanAgg = accumarray(aggindex, janTairMean, [numel(sites) 1], @nanmean);
janTs20meanAgg = accumarray(aggindex, janTs20mean, [numel(sites) 1], @nanmean);

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
latAgg = accumarray(aggindex, lat, [numel(sites) 1], @mean);
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
selectState = 'UT';
selectIDs = siteIDs(strcmpi(states, selectState));
st_test = ismember(sites_cl, selectIDs);

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
h = figure(fignum);
set(h, 'Name','(UT-Agg) January Air/Soil T gradients',...
    'DefaultAxesFontSize',18, 'DefaultTextFontSize', 16);

% Aggregate the mean, sd, AND the sites_cl list
janTairSdAgg = accumarray(aggindex, janTairSd, [numel(sites) 1], @mean);
janTs20sdAgg = accumarray(aggindex, janTs20sd, [numel(sites) 1], @mean);
sites_clAgg = accumarray(aggindex, sites_cl, [numel(sites) 1], @mean);

% Get selected state site out of aggregated data
st_testAgg = ismember(sites_clAgg, selectIDs);

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
    'MarkerEdgeColor', 'k','MarkerSize', 8);
hold on;
xfit = linspace(1692, 3400);
[b,bint,resid,rint,stats] = shregress(y1, [x ones(size(x))]);
handles(2) = plot(xfit, polyval(b, xfit), '--k');
if stats(3) < 0.01
    text(1750, -7, ['r^2 = ' num2str(stats(1),2) ', p < 0.01']);
else
    text(1750, -7, ['r^2 = ' num2str(stats(1),2)...
        ', p = ' num2str(stats(3),2)]);
end

% Plot Tair example error bar in an inset box
handles(3) = errorbar(1625, -4.7, meanStdDev,'o',...
    'MarkerFaceColor',[0.7 0.7 0.7], 'Color', [0.7 0.7 0.7],...
    'MarkerSize', 8);
% Create rectangle
annotation(h, 'rectangle', [0.154 0.146 0.05 0.55], 'FaceColor','flat');

% Plot Jan mean Tsoil and regression
handles(4) = errorbar(x, y2, sd2, 'ok', 'MarkerFaceColor', 'w',...
    'Color', [0.5 0.5 0.5], 'MarkerEdgeColor', 'k','MarkerSize', 8);
xfit = linspace(1600, 3400);
[b,bint,resid,rint,stats] = shregress(y2, [x ones(size(x))]);
handles(5) = plot(xfit, polyval(b, xfit), '--k');
if stats(3) < 0.01
    text(1750, 3, ['r^2 = ' num2str(stats(1),2) ', p < 0.01']);
else
    text(1750, 3, ['r^2 = ' num2str(stats(1),2)...
        ', p = ' num2str(stats(3),2)]);
end

xlim([1500 3500]);
xlabel('Elevation(m)');
ylabel('^oC');
legend(handles([1 4]), {'T_{air}',...
    'T_{soil} (20cm)'});
title('Utah - January Tair/Tsoil (aggregated)');
% Plot moist adiabatic lapse rate
% plot([900 3500], [4, -9], ':k')

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

% Aggregate the mean and sd
julTs20sdAgg = accumarray(aggindex, julTs20sd, [numel(sites) 1], @mean);

% Assign x and y variables
x = elevAgg(st_testAgg);
y1 = janTs20meanAgg(st_testAgg);
sd1 = janTs20sdAgg(st_testAgg);
y2 = julTs20meanAgg(st_testAgg);
sd2 = julTs20sdAgg(st_testAgg);

fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name','(UT-Agg) January/July Soil T gradients',...
    'DefaultAxesFontSize',18, 'DefaultTextFontSize', 16);
% Plot January Ts by elevation and regression
handles(1) = errorbar(x, y1, sd1, 'ok', 'MarkerFaceColor', 'w',...
    'MarkerSize', 8);
hold on;
xfit = linspace(1600, 3400);
[b,bint,resid,rint,stats] = shregress(y1, [x ones(size(x))]);
handles(2) = plot(xfit, polyval(b, xfit), '--k');
if stats(3) < 0.01
    text(1650, 3, ['r^2 = ' num2str(stats(1),2) ', p < 0.01']);
else
    text(1650, 3, ['r^2 = ' num2str(stats(1),2)...
        ', p = ' num2str(stats(3),2)]);
end
% Plot July Ts by elevation and regression
handles(3) = errorbar(x, y2, sd2, 'ok', 'MarkerFaceColor', 'k',...
    'MarkerSize', 8);
[b,bint,resid,rint,stats] = shregress(y2, [x ones(size(x))]);
handles(4) = plot(xfit, polyval(b, xfit), '--k');
if stats(3) < 0.01
    text(1650, 18.5, ['r^2 = ' num2str(stats(1),2) ', p < 0.01']);
else
    text(1650, 18.5, ['r^2 = ' num2str(stats(1),2)...
        ', p = ' num2str(stats(3),2)]);
end
% Plot moist adiabatic lapse rate
%plot([1600 3400], [18, 9], ':k');
xlim([1500 3500]);
legend(handles([1 3]), {'January', 'July'});
xlabel('Elevation (m)');
ylabel('20cm T_{soil} (^oC)')

junk = 99;