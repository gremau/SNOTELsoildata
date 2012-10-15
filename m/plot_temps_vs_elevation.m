% plot_t_elevgradients.m
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
%addpath('../m/');
addpath('/home/greg/data/code_resources/m_common/nanstuff/');

% Set data path and file name, read in file
rawdatapath = '../rawdata/';
processeddatapath = '../processed_data/';

% Load lists of sites with data in the daily/hourly data directory
dailysites = sortrows(csvread('../rawdata/allsensors_daily/filelist.txt'));
soilsites = sortrows(csvread('../rawdata/soilsensors_hourly/filelist.txt'));
sites = unique(soilsites(:,1));

% LOAD the data (can switch between daily/hourly data here)
climData = csvread([processeddatapath 'wyear_climatesummary.txt']);
tsData = csvread([processeddatapath 'wyear_soiltempsummary_hourly.txt']);
% tsData = csvread([processeddatapath 'wyear_soiltempsummary_daily.txt']);

% Get a subset of climData that corresponds with available soildata
[matchtest, idx] = ismember(climData(:, 1:2), tsData(:, 1:2), 'rows');
soilClim = climData(matchtest, :);
% matchtest2 = ismember(tsData(:, 1:2), soilClim(:, 1:2), 'rows');

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
maxswe = soilClim(:,3)*25.4;
maat = soilClim(:, 74); % Mean wateryear air temp
onsetdoy = soilClim(:, 6);
meltdoy = soilClim(:, 7);
totaldaysSC = soilClim(:, 9);
octAirTmean = soilClim(:,50);
octAirTsd = soilClim(:, 51);
novAirTmean = soilClim(:, 52);
novAirTsd = soilClim(:, 53);
decAirTmean = soilClim(:, 54);
decAirTsd = soilClim(:, 55);
janAirTmean = soilClim(:, 56);
janAirTsd = soilClim(:, 57);
febAirTmean = soilClim(:, 58);
febAirTsd = soilClim(:, 59);
marAirTmean = soilClim(:, 60);
marAirTsd = soilClim(:, 61);
aprAirTmean = soilClim(:, 62);
aprAirTsd = soilClim(:, 63);
mayAirTmean = soilClim(:, 64);
mayAirTsd = soilClim(:, 65);
junAirTmean = soilClim(:, 66);
junAirTsd = soilClim(:, 67);
julAirTmean = soilClim(:, 68);
julAirTsd = soilClim(:, 69);
augAirTmean = soilClim(:, 70);
augAirTsd = soilClim(:, 71);
sepAirTmean = soilClim(:, 72);
sepAirTsd = soilClim(:, 73);

elev = soilClim(:, 82);
lat = soilClim(:, 83);
lon = soilClim(:, 84);
ltMeanSWE = soilClim(:, 85);
ltMeanPrecip = soilClim(:, 86);

oct5cmSTmean = tsData(:, 3);
oct5cmSTsd = tsData(:, 4);
oct20cmSTmean = tsData(:, 5);
oct20cmSTsd = tsData(:, 6);
oct50cmSTmean = tsData(:, 7);
oct50cmSTsd = tsData(:, 8);
nov5cmSTmean = tsData(:, 9);
nov5cmSTsd = tsData(:, 10);
nov20cmSTmean = tsData(:, 11);
nov20cmSTsd = tsData(:, 12);
nov50cmSTmean = tsData(:, 13);
nov50cmSTsd = tsData(:, 14);
dec5cmSTmean = tsData(:, 15);
dec5cmSTsd = tsData(:, 16);
dec20cmSTmean = tsData(:, 17);
dec20cmSTsd = tsData(:, 18);
dec50cmSTmean = tsData(:, 19);
dec50cmSTsd = tsData(:, 20);
jan5cmSTmean = tsData(:, 21);
jan5cmSTsd = tsData(:, 22);
jan20cmSTmean = tsData(:, 23);
jan20cmSTsd = tsData(:, 24);
jan50cmSTmean = tsData(:, 25);
jan50cmSTsd = tsData(:, 26);
feb5cmSTmean = tsData(:, 27);
feb5cmSTsd = tsData(:, 28);
feb20cmSTmean = tsData(:, 29);
feb20cmSTsd = tsData(:, 30);
feb50cmSTmean = tsData(:, 31);
feb50cmSTsd = tsData(:, 32);
mar5cmSTmean = tsData(:, 33);
mar5cmSTsd = tsData(:, 34);
mar20cmSTmean = tsData(:, 35);
mar20cmSTsd = tsData(:, 36);
mar50cmSTmean = tsData(:, 37);
mar50cmSTsd = tsData(:, 38);
apr5cmSTmean = tsData(:, 39);
apr5cmSTsd = tsData(:, 40);
apr20cmSTmean = tsData(:, 41);
apr20cmSTsd = tsData(:, 42);
apr50cmSTmean = tsData(:, 43);
apr50cmSTsd = tsData(:, 44);
may5cmSTmean = tsData(:, 45);
may5cmSTsd = tsData(:, 46);
may20cmSTmean = tsData(:, 47);
may20cmSTsd = tsData(:, 48);
may50cmSTmean = tsData(:, 49);
may50cmSTsd = tsData(:, 50);
jun5cmSTmean = tsData(:, 51);
jun5cmSTsd = tsData(:, 52);
jun20cmSTmean = tsData(:, 53);
jun20cmSTsd = tsData(:, 54);
jun50cmSTmean = tsData(:, 55);
jun50cmSTsd = tsData(:, 56);
jul5cmSTmean = tsData(:, 57);
jul5cmSTsd = tsData(:, 58);
jul20cmSTmean = tsData(:, 59);
jul20cmSTsd = tsData(:, 60);
jul50cmSTmean = tsData(:, 61);
jul50cmSTsd = tsData(:, 62);
aug5cmSTmean = tsData(:, 63);
aug5cmSTsd = tsData(:, 64);
aug20cmSTmean = tsData(:, 65);
aug20cmSTsd = tsData(:, 66);
aug50cmSTmean = tsData(:, 67);
aug50cmSTsd = tsData(:, 68);
sep5cmSTmean = tsData(:, 69);
sep5cmSTsd = tsData(:, 70);
sep20cmSTmean = tsData(:, 71);
sep20cmSTsd = tsData(:, 72);
sep50cmSTmean = tsData(:, 73);
sep50cmSTsd = tsData(:, 74);

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

% PLOTS
%----------------------------------
% FIG 1 - January and July Tair/Tsoil gradients - All SNOTEL/years
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name','Jan/July Air and soil temperature gradients, All SNOTEL');

subplot (2, 2, 1)
plot(elev, julAirTmean, 'ok', 'MarkerFaceColor', 'r');
hold on
plot(elev, jul20cmSTmean, 'ok', 'MarkerFaceColor', ...
    [0.7 0.7 0.7]);
xlabel('Elevation(m)');
ylabel('Mean July temp (^oC)');
title('July');
legend('Air', 'Soil(20cm)');

subplot (2, 2, 2)
plot(elev, janAirTmean, 'ok', 'MarkerFaceColor', 'b');
hold on
plot(elev, jan20cmSTmean, 'ok', 'MarkerFaceColor', ...
    [0.7 0.7 0.7]);
% Plot moist adiabatic lapse rate
plot([900 3500], [4, -9], ':k')
xlabel('Elevation(m)');
ylabel('Mean temp (Celsius)');
title('January');
legend('Air', 'Soil', 'Moist adiabatic lapse( 5^oC/km)');

% Aggregate by wateryear using accumarray
elevAgg = accumarray(aggindex, elev, [numel(sites) 1], @mean);
julAirTmeanAgg = accumarray(aggindex, julAirTmean, [numel(sites) 1], @nanmean);
jul20cmSTmeanAgg = accumarray(aggindex, jul20cmSTmean, [numel(sites) 1], @nanmean);
janAirTmeanAgg = accumarray(aggindex, janAirTmean, [numel(sites) 1], @nanmean);
jan20cmSTmeanAgg = accumarray(aggindex, jan20cmSTmean, [numel(sites) 1], @nanmean);

subplot (2, 2, 3)
plot(elevAgg, julAirTmeanAgg, 'ok', 'MarkerFaceColor', 'r');
hold on
plot(elevAgg, jul20cmSTmeanAgg, 'ok', 'MarkerFaceColor', ...
    [0.7 0.7 0.7]);
xlabel('Elevation(m)');
ylabel('Mean July temp (^oC)');
title('July (mean of all wateryears)');
legend('Air', 'Soil(20cm)');

subplot (2, 2, 4)
plot(elevAgg, janAirTmeanAgg, 'ok', 'MarkerFaceColor', 'b');
hold on
plot(elevAgg, jan20cmSTmeanAgg, 'ok', 'MarkerFaceColor', ...
    [0.7 0.7 0.7]);
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
errorbar(elev, snowcovMeanST20cm, snowcovStdST20cm, 'b.');
hold on
errorbar(elev, snowfreeMeanST20cm, snowfreeStdST20cm, 'r.');
title('Snowcovered and snowfree temperature');
legend('Snow-covered', 'Snow-free');
ylabel('20cm soil temp');
xlabel('Elevation (m)');

subplot(2,2,2);
errorbar(maxswe, snowcovMeanST20cm, snowcovStdST20cm, 'b.');
hold on
errorbar(maxswe, snowfreeMeanST20cm, snowfreeStdST20cm, 'r.');
title('Snowcovered and snowfree temperature');
legend('Snow-covered', 'Snow-free');
ylabel('20cm Ts (^oC)');
xlabel('Peak SWE (mm)');

% Plot January and July soil temps by elevation
subplot(2,2,3);
errorbar(elev, jan20cmSTmean,jan20cmSTsd, 'ko', 'MarkerFaceColor', 'w');
hold on;
errorbar(elev, jul20cmSTmean, jul20cmSTsd, 'ko', 'MarkerFaceColor', 'k');
% Plot moist adiabatic lapse rate
plot([1600 3400], [18, 9], ':k');
legend('January', 'July', 'Moist adiabatic lapse (5^oC/km)');
xlabel('Elevation (m)');
ylabel('20cm Ts (^oC)')

% Plot January and July soil temps by elevation
subplot(2,2,4);
errorbar(maxswe, jan20cmSTmean,jan20cmSTsd, 'ko', 'MarkerFaceColor', 'w');
hold on;
errorbar(maxswe, jul20cmSTmean, jul20cmSTsd, 'ko', 'MarkerFaceColor', 'k');
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
plot(elevAgg, janAirTmeanAgg, '.b');
plot(elevAgg, julAirTmeanAgg, '.m');
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
plot(latAgg, julAirTmeanAgg, '.', 'Color', [0.4 0.4 0.4]);
plot(latAgg, janAirTmeanAgg, '.k');
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
    eval(['y1 = ' months(i,:) 'AirTmean;']);
    eval(['y2 = ' months(i,:) '5cmSTmean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
%     [~, rsq, xfit, yfit] = fitline(x, y1, polyorder, [0, 1500]);
%     plot(xfit, yfit, '-k', 'LineWidth', 2);
%     text(700, 8, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
    plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    title(months(i,:));
    xlim([800, 4000]); ylim([-15, 30]);
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('^oC')
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
    
    subplot(4,3, i+3)
    eval(['y1 = ' months(i+3,:) 'AirTmean;']);
    eval(['y2 = ' months(i+3,:) '20cmSTmean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    title(months(i+3,:));
    xlim([800, 4000]); ylim([-15, 30]);
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('^oC')
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
    
    subplot(4,3, i+6);
    eval(['y1 = ' months(i+6,:) 'AirTmean;']);
    eval(['y2 = ' months(i+6,:) '50cmSTmean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    title(months(i+6,:));
    xlim([800, 4000]); ylim([-15, 30]);
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('^oC')
    elseif i>1
        set(gca, 'YTickLabel', '');
    end

    subplot(4,3, i+9);
    eval(['y1 = ' months(i+9,:) 'AirTmean;']);
    eval(['y2 = ' months(i+9,:) '50cmSTmean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    title(months(i+9,:));
    xlim([800, 4000]); ylim([-15, 30]);
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
    eval(['y1 = ' months(i,:) 'AirTmean;']);
    eval(['y2 = ' months(i,:) '5cmSTmean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
%     [~, rsq, xfit, yfit] = fitline(x, y1, polyorder, [0, 1500]);
%     plot(xfit, yfit, '-k', 'LineWidth', 2);
%     text(700, 8, ['r^2 = ' num2str(rsq, 2)]); % r^2 values
    plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    xlim([800, 4000]); ylim([-10, 10]);
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('^oC')
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
    
    subplot(3, 5, i+5)
    eval(['y2 = ' months(i,:) '20cmSTmean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    xlim([800, 4000]); ylim([-10, 10]);
    set(gca, 'XTickLabel', '');
    if i==1
        ylabel('^oC')
    elseif i>1
        set(gca, 'YTickLabel', '');
    end
    
    subplot(3, 5, i+10);
    eval(['y2 = ' months(i,:) '50cmSTmean;']);
    plot(x, y1, '.', 'Color', [0.7,0.7,0.7]);
    hold on;
    plot(x, y2, '.', 'Color', [0.4,0.4,0.4]);
    xlim([800, 4000]); ylim([-10, 10]);
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
plot(ltMeanSWE(swetest), julAirTmean(swetest), '.k');
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
plot(elev(temptest), julAirTmean(temptest), '.k');
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
plot(elev, decAirTmean, 'ok', 'MarkerFaceColor', ...
    [0.7 0.7 0.7]);
hold on
plot(elev, dec20cmSTmean, 'ok', 'MarkerFaceColor', 'r');
% Plot moist adiabatic lapse rate
plot([900 3500], [4, -9], ':k')
xlabel('Elevation(m)');
ylabel('Mean temp (Celsius)');
legend('Air', 'Soil(20cm)', 'Moist adiabatic lapse( 5^oC/km)');

subplot (1,2,2);
plot(elev, decAirTmean, 'ok', 'MarkerFaceColor', ...
    [0.7 0.7 0.7]);
hold on
plot(elev, dec5cmSTmean, 'ok', 'MarkerFaceColor', 'r');
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
errorbar(elev, decAirTmean, decAirTsd, ...
    'ok', 'MarkerFaceColor', [0.7 0.7 0.7]);
hold on
errorbar(elev, dec20cmSTmean, dec20cmSTsd, ...
    'ok', 'MarkerFaceColor', 'r');
% Plot moist adiabatic lapse rate
plot([900 3500], [4, -9], ':k')
xlabel('Elevation(m)');
ylabel('Mean temp (^oC)');
legend('Air', 'Soil (20cm)', 'Moist adiabatic lapse( 5^oC/km)');
title('Air & soil temp.');

subplot(1,2,2);
offset = abs(decAirTmean - dec20cmSTmean);
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

selstate = 'UT';
selsites = siteIDs(strcmpi(states, 'UT'));
st_test = ismember(soilClim(:,1), selsites);
% statetest = 

subplot (1,1,1);
errorbar(elev(st_test), janAirTmean(st_test), janAirTsd(st_test), ...
    'ok', 'MarkerFaceColor', [0.7 0.7 0.7], 'Color', [0.7 0.7 0.7], ...
    'MarkerEdgeColor', 'k');
hold on
errorbar(elev(st_test), jan20cmSTmean(st_test), jan20cmSTsd(st_test), ...
    'ok', 'MarkerFaceColor', 'r', 'Color', [0.7 0.7 0.7], ...
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
set(h, 'Name','(UT-Agg) January Air/Soil T gradients');

st_test = ismember(unique(soilClim(:,1)), selsites);
% Aggregate the mean and sd
janAirTsdAgg = accumarray(aggindex, janAirTsd, [numel(sites) 1], @mean);
jan20cmSTsdAgg = accumarray(aggindex, jan20cmSTsd, [numel(sites) 1], @mean);

subplot (1,1,1);
errorbar(elevAgg(st_test), janAirTmeanAgg(st_test), janAirTsdAgg(st_test), ...
    'ok', 'MarkerFaceColor', [0.7 0.7 0.7], 'Color', [0.7 0.7 0.7], ...
    'MarkerEdgeColor', 'k');
hold on
errorbar(elevAgg(st_test), jan20cmSTmeanAgg(st_test), jan20cmSTsdAgg(st_test), ...
    'ok', 'MarkerFaceColor', 'r', 'Color', [0.7 0.7 0.7], ...
    'MarkerEdgeColor', 'k');
% Plot moist adiabatic lapse rate
plot([900 3500], [4, -9], ':k')
xlim([1500 3500]);
xlabel('Elevation(m)');
ylabel('Mean temp (^oC)');
legend('Air', 'Soil (20cm)', 'Moist adiabatic lapse( 5^oC/km)');
title('Utah - January Tair/Tsoil (aggregated)');


% -------------------------------------------------------------
% FIG 11 - January/July soil T gradients - Utah data only

fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name','(UT) January/July Soil T gradients');
% Plot January and July soil temps by elevation
errorbar(elev(st_test), jan20cmSTmean(st_test),...
    jan20cmSTsd(st_test), 'ko', 'MarkerFaceColor', 'w');
hold on;
errorbar(elev(st_test), jul20cmSTmean(st_test),...
    jul20cmSTsd(st_test), 'ko', 'MarkerFaceColor', 'k');
% Plot moist adiabatic lapse rate
plot([1600 3400], [18, 9], ':k');
legend('January', 'July', 'Moist adiabatic lapse (5^oC/km)');
xlabel('Elevation (m)');
ylabel('20cm Ts (^oC)')

% -------------------------------------------------------------
% FIG 12 - January/July soil T gradients - AggregatedUtah data only

% Aggregate the mean and sd
%julAirTsdAgg = accumarray(aggindex, julAirTsd, [numel(sites) 1], @mean);
jul20cmSTsdAgg = accumarray(aggindex, jul20cmSTsd, [numel(sites) 1], @mean);

fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name','(UT-Agg) January/July Soil T gradients');
% Plot January and July soil temps by elevation
errorbar(elevAgg(st_test), jan20cmSTmeanAgg(st_test),...
    jan20cmSTsdAgg(st_test), 'ko', 'MarkerFaceColor', 'w');
hold on;
errorbar(elevAgg(st_test), jul20cmSTmeanAgg(st_test),...
    jul20cmSTsdAgg(st_test), 'ko', 'MarkerFaceColor', 'k');
% Plot moist adiabatic lapse rate
plot([1600 3400], [18, 9], ':k');
legend('January', 'July', 'Moist adiabatic lapse (5^oC/km)');
xlabel('Elevation (m)');
ylabel('20cm Ts (^oC)')

junk = 99;
