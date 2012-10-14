% plot_temps_vs_elevation.m
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
addpath('/home/greg/data/programming_resources/m_common/');

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

%-----------------------------
% Assign variables
decAirTmean = soilClim(:, 54);
decAirTsd = soilClim(:, 55);
julAirTmean = soilClim(:, 68);
julAirTsd = soilClim(:, 69);

elev = soilClim(:, 82);
lat = soilClim(:, 83);
lon = soilClim(:, 84);
ltMeanSWE = soilClim(:, 85);
ltMeanPrecip = soilClim(:, 86);
dec5cmSTmean = tsData(:, 15);
dec5cmSTsd = tsData(:, 16);
dec20cmSTmean = tsData(:, 17);
dec20cmSTsd = tsData(:, 18);
dec50cmSTmean = tsData(:, 19);
dec50cmSTsd = tsData(:, 20);
jul5cmSTmean = tsData(:, 57);
jul5cmSTsd = tsData(:, 58);
jul20cmSTmean = tsData(:, 59);
jul20cmSTsd = tsData(:, 60);
jul50cmSTmean = tsData(:, 61);
jul50cmSTsd = tsData(:, 62);


% Ask user for month number and state
monthsel = str2double(input('Which month (1-12)?: ', 's'));
statesel = input(...
    'Which state ("AZ, CO, ID, MT, NM, NV, UT, WY, or all")?: ', 's');
monthlabels = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sept' 'Oct'...
    'Nov' 'Dec'};
monthlabel = monthlabels{monthsel};

%Load list of sites in the daily data directory
havedata = csvread([rawdatapath 'allsensors_daily/filelist.txt']);
havedata = unique(havedata(:, 1));

% Load 30 year average data
avgSWE = load7100Avg('swe');
avgPrecip = load7100Avg('precip');

% Get an inventory of  sites and their elevations from inventory file
% Create format string (station,elev,cdbs_id only here)
formatstr = '%*s%*s%*s%*s%s%f%f%f%f';
fid = fopen([processeddatapath 'SNOTELinventory.csv']);
inventorycell = textscan(fid, formatstr,'Headerlines', 1, 'Delimiter', ',');
fclose(fid);

states = inventorycell{1};
sitesarray = [inventorycell{2}, inventorycell{5}*0.3048];
clear inventorycell;

% Remove site rows that are not part of the selected state
if strcmpi(statesel, 'all')
    sitesarray = sitesarray;
else
    stateseltest = strcmpi(states, statesel);
    sitesarray = sitesarray(stateseltest, :);
end

% Filter out site rows for which there is no data in the data folder
havedatatest = ismember(sitesarray(:, 1), havedata);
sitesarray = sitesarray(havedatatest, :);

% Load data for each site into a large cellarray
datacell = cell(length(sitesarray), 1);
for i = 1:length(sitesarray);
    m = loadsnotel(sitesarray(i,1), 'daily', 'exclude');
    dailysite_datevec = datevec(m{2}, 'yyyy-mm-dd');
    % Create logical test for July and selected month
    julytest = dailysite_datevec(:,2)==7;
    monthtest = dailysite_datevec(:,2)==monthsel;
    % Append mean July and selected month airtemps on to sitesarray
    sitesarray(i, 3) = nanmean(m{7}(julytest));% July mean max air temp
    sitesarray(i, 4) = nanmean(m{9}(julytest));% july mean avg air temp
    sitesarray(i, 5) = nanmean(m{15}(julytest));% july mean 20cm soil temp
    sitesarray(i, 6) = nanmean(m{7}(monthtest));% mean max air temp
    sitesarray(i, 7) = nanmean(m{9}(monthtest));% mean avg air temp
    sitesarray(i, 8) = nanmean(m{15}(monthtest));% mean 20cm soil temp
    sitesarray(i, 9) = nanmean(m{14}(monthtest));% mean 5cm soil temp
    %sitesarray(i, 12) = nanstd(m{9}(monthtest));% stdev avg air temp
    %sitesarray(i, 13) = nanstd(m{15}(monthtest));% stdev 20cm soil temp
    % Test if there is an average SWE or Precip
    avgSWEtest = avgSWE(:,1)==sitesarray(i,1);
    avgPreciptest = avgPrecip(:,1)==sitesarray(i,1);
%     twoavg = sum(avgSWEtest) + sum(avgPreciptest);
    if sum(avgSWEtest)==1
        sitesarray(i, 11) = avgSWE(avgSWEtest, 27); % AVG peakswe
    else
        sitesarray(i, 11) = nan;
    end
    if sum(avgPreciptest)==1
        sitesarray(i, 10) = avgPrecip(avgPreciptest, 15); % AVG precip
    else
        sitesarray(i, 10) = nan; %AVG precip
    end
    sitesarray(i, 12) = nanstd(m{9}(monthtest));% stdev avg air temp
    sitesarray(i, 13) = nanstd(m{15}(monthtest));% stdev 20cm soil temp
end

% PLOTS
%----------------------------------
% FIG 1 - Full  gradients
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name',['Mean air and soil temps by elevation, ' statesel ...
    ' SNOTEL sites']);

subplot (2, 2, 1)
plot(sitesarray(:,2), sitesarray(:,4), '.m');
hold on
plot(sitesarray(:,2), sitesarray(:,5), '.b');
xlabel('Elevation(m)');
ylabel('Mean July temp (^oC)');
title('JULY - Air and soil temp vs elevation');
legend('Air', 'Soil(20cm, 12am)');

subplot (2, 2, 2)
plot(sitesarray(:,2), sitesarray(:,7), 'ok', 'MarkerFaceColor', ...
    [0.7 0.7 0.7]);
hold on
plot(sitesarray(:,2), sitesarray(:,8), 'ok', 'MarkerFaceColor', 'r');
% Plot moist adiabatic lapse rate
plot([900 3500], [4, -9], ':k')
xlabel('Elevation(m)');
ylabel('Mean temp (Celsius)');
title([ monthlabel ' - Air and soil temp vs elevation']);
legend('Air', 'Soil(20cm)', 'Moist adiabatic lapse( 5^oC/km)');

subplot (2, 2, 3)
plot(sitesarray(:,2), sitesarray(:,3), '.r');
hold on
plot(sitesarray(:,2), sitesarray(:,6), '.b');
xlabel('Elevation(m)');
ylabel('MAX air temp (Celsius)');
title('Max air temp vs elevation');
legend('July', [ monthlabel, '(20cm, 12am)']);

subplot (2, 2, 4)
plot(sitesarray(:,2), sitesarray(:,10), '.k');
hold on
plot(sitesarray(:,2), sitesarray(:,11), '.b');
xlabel('Elevation(m)');
ylabel('mm of water (Celsius)');
title('Precip and SWE');
legend('Total precip', 'Peak SWE');

%----------------------------------
% FIG 1 - Full  gradients
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name',['Mean air and soil temps by elevation, ' statesel ...
    ' SNOTEL sites']);

subplot (2, 2, 1)
plot(elev, julAirTmean, '.m');
hold on
plot(elev, jul20cmSTmean, '.b');
xlabel('Elevation(m)');
ylabel('Mean July temp (^oC)');
title('JULY - Air and soil temp vs elevation');
legend('Air', 'Soil(20cm, 12am)');

subplot (2, 2, 2)
plot(elev, decAirTmean, 'ok', 'MarkerFaceColor', ...
    [0.7 0.7 0.7]);
hold on
plot(elev, dec20cmSTmean, 'ok', 'MarkerFaceColor', 'r');
% Plot moist adiabatic lapse rate
plot([900 3500], [4, -9], ':k')
xlabel('Elevation(m)');
ylabel('Mean temp (Celsius)');
title([ monthlabel ' - Air and soil temp vs elevation']);
legend('Air', 'Soil(20cm)', 'Moist adiabatic lapse( 5^oC/km)');

subplot (2, 2, 3)
plot(sitesarray(:,3), sitesarray(:,3), '.r');
hold on
plot(sitesarray(:,3), sitesarray(:,6), '.b');
xlabel('Elevation(m)');
ylabel('MAX air temp (Celsius)');
title('Max air temp vs elevation');
legend('July', [ monthlabel, '(20cm, 12am)']);

subplot (2, 2, 4)
plot(elev, ltMeanPrecip, '.k');
hold on
plot(elev, ltMeanSWE, '.b');
xlabel('Elevation(m)');
ylabel('mm of water (Celsius)');
title('Precip and SWE');
legend('Total precip', 'Peak SWE');


%----------------------------------
% FIG 1 - Full  gradients - aggregate all years
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name',['Mean air and soil temps by elevation, ' statesel ...
    ' SNOTEL sites']);

[vals, ~, valindex] = unique(soilClim(:,1));
aggindex = [valindex ones(size(soilClim(:,1)))];

subplot (2, 2, 1)
% First accumilate the data with accumarray
elevAgg = accumarray(aggindex, elev, [numel(sites) 1], @mean);
julAirTmeanAgg = accumarray(aggindex, julAirTmean, [numel(sites) 1], @nanmean);
jul20cmSTmeanAgg = accumarray(aggindex, jul20cmSTmean, [numel(sites) 1], @nanmean);

plot(elevAgg, julAirTmeanAgg, '.m');
hold on
plot(elevAgg, jul20cmSTmeanAgg, '.b');
xlabel('Elevation(m)');
ylabel('Mean July temp (^oC)');
title('JULY - Air and soil temp vs elevation');
legend('Air', 'Soil(20cm, 12am)');

subplot (2, 2, 2)

decAirTmeanAgg = accumarray(aggindex, decAirTmean, [numel(sites) 1], @nanmean);
dec20cmSTmeanAgg = accumarray(aggindex, dec20cmSTmean, [numel(sites) 1], @nanmean);

plot(elevAgg, decAirTmeanAgg, 'ok', 'MarkerFaceColor', ...
    [0.7 0.7 0.7]);
hold on
plot(elevAgg, dec20cmSTmeanAgg, 'ok', 'MarkerFaceColor', 'r');
% Plot moist adiabatic lapse rate
plot([900 3500], [4, -9], ':k')
xlabel('Elevation(m)');
ylabel('Mean temp (Celsius)');
title([ monthlabel ' - Air and soil temp vs elevation']);
legend('Air', 'Soil(20cm)', 'Moist adiabatic lapse( 5^oC/km)');

subplot (2, 2, 3)
plot(sitesarray(:,3), sitesarray(:,3), '.r');
hold on
plot(sitesarray(:,3), sitesarray(:,6), '.b');
xlabel('Elevation(m)');
ylabel('MAX air temp (Celsius)');
title('Max air temp vs elevation');
legend('July', [ monthlabel, '(20cm, 12am)']);

subplot (2, 2, 4)
plot(elev, ltMeanPrecip, '.k');
hold on
plot(elev, ltMeanSWE, '.b');
xlabel('Elevation(m)');
ylabel('mm of water (Celsius)');
title('Precip and SWE');
legend('Total precip', 'Peak SWE');

%--------------------------------------------------------
% FIG 2 - Restricted gradients
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name','Restricted gradients');

%Select peak SWE gradient
swetest = sitesarray(:, 11)<600 & sitesarray(:, 11)>400;
swesubset = sitesarray(swetest, :);

subplot (2, 2, 1)
plot(swesubset(:,11), swesubset(:,4), '.k');
xlabel('Mean peak SWE(mm)');
ylabel('Mean July temp (Celsius)');
title('Constant SWE - July air temp');
%legend('Air', 'Soil(20cm, 12am)');

subplot (2, 2, 2)
plot(swesubset(:,11), swesubset(:,2), '.k');
xlabel('Mean peak SWE(mm)');
ylabel('elevation(m)');
title('Constant SWE - elevation');

temptest = sitesarray(:, 2)<2550 & sitesarray(:, 2)>2350;
tempsubset = sitesarray(temptest, :);

subplot (2, 2, 3)
plot(tempsubset(:,2), tempsubset(:,4), '.k');
xlabel('elevation(m)');
ylabel('Mean July temperature (Celsius');
title('Constant elevation - July temperature');

subplot (2, 2, 4)
plot(tempsubset(:,2), tempsubset(:,11), '.k');
xlabel('elevation(m)');
ylabel('Mean peak SWE(mm)');
title('Constant elevation - Peak SWE');


% ------------------------------------------------------
% FIG 3 - Air temp gradient
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name', [ monthlabel ' air temp vs elevation']);

plot(sitesarray(:,2), sitesarray(:,7), '.m');
hold on
% Plot moist adiabatic lapse rate
plot([1200 3500], [-1, -8.2], ':k')
% plot(sitesarray(:,2), sitesarray(:,7), '.b');
xlabel('elevation(m)');
ylabel('Mean air temp (Celsius)');
title([monthlabel ' - Air temp vs elevation']);

%----------------------------------------------------------------------
% FIG 4
fignum = fignum+1;    
h = figure(fignum);

subplot (1,2,1);
plot(sitesarray(:,2), sitesarray(:,7), 'ok', 'MarkerFaceColor', ...
    [0.7 0.7 0.7]);
hold on
plot(sitesarray(:,2), sitesarray(:,8), 'ok', 'MarkerFaceColor', 'r');
% Plot moist adiabatic lapse rate
plot([900 3500], [4, -9], ':k')
xlabel('Elevation(m)');
ylabel('Mean temp (Celsius)');
legend('Air', 'Soil(20cm)', 'Moist adiabatic lapse( 5^oC/km)');

subplot (1,2,2);
plot(sitesarray(:,2), sitesarray(:,7), 'ok', 'MarkerFaceColor', ...
    [0.7 0.7 0.7]);
hold on
plot(sitesarray(:,2), sitesarray(:,9), 'ok', 'MarkerFaceColor', 'r');
% Plot moist adiabatic lapse rate
plot([900 3500], [4, -9], ':k')
xlabel('Elevation(m)');
ylabel('Mean temp (Celsius)');
legend('Air', 'Soil(5cm)', 'Moist adiabatic lapse( 5^oC/km)');



% -------------------------------------------------------------
% FIG 5
fignum = fignum+1;    
h = figure(fignum);

subplot (1,2,1);
errorbar(sitesarray(:,2), sitesarray(:,7), sitesarray(:,12), ...
    'ok', 'MarkerFaceColor', [0.7 0.7 0.7]);
hold on
errorbar(sitesarray(:,2), sitesarray(:,8), sitesarray(:,13), ...
    'ok', 'MarkerFaceColor', 'r');
% Plot moist adiabatic lapse rate
plot([900 3500], [4, -9], ':k')
xlabel('Elevation(m)');
ylabel('Mean temp (^oC)');
legend('Air', 'Soil (20cm)', 'Moist adiabatic lapse( 5^oC/km)');
title([ monthlabel ' air & soil temp.']);

subplot(1,2,2);
offset = abs(sitesarray(:,7) - sitesarray(:,8));
plot(sitesarray(:,2), offset(:,1), 'ok', 'MarkerFaceColor', 'b');
xlabel('Elevation (m)');
ylabel('^oC');
title('Air-Soil temp. offset');
%legend('Air-Soil Offset');

junk = 99;
