% plot_temps_vs_elevation.m
%
% version 3: GM, 120216
% Made it a user-selectable month (still need to change variable names)


clear;          % clear memory
close all;      % clear any figures
fignum=0;       % used to increment figure number for plots
%addpath('../m/');
addpath('/home/greg/data/programming_resources/m_common/');

% Set data path and file name, read in file
datapath = '../rawdata/';

% Ask user for month number
month = str2double(input('Which month (1-12)?: ', 's'));
monthLabels = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sept' 'Oct'...
    'Nov' 'Dec'};

%Load list of sites in the daily data directory
havedata = unique(dlmread([datapath 'allsensors_daily/sitelist.txt']));

% Load list of sites with data in the daily data directory
dailyDataSites = sortrows(csvread([datapath 'allsensors_daily/sitelist.txt']));

% Load list of bad data years for all sites
badDataYears = sortrows(csvread([datapath 'allsensors_daily/baddata.txt'], 1, 0));

% Load 30 year average data
averages = dlmread([datapath 'longterm_averages/'...
    '7100_avgprecipswe_snotel.csv'], ',', 1, 4);

% Generate list of sites and their elevations from inventory file
% Create format string (station,elev,cdbs_id only here)
formatstr = '%*s%f%*s%*s%*s%*s%f%*s%s%*s%*s%*s%*s%*s%*s%*u';
fid = fopen([datapath 'station_inventory/UT_soilstations.csv']);
%fid = fopen([datapath 'station_inventory/merged_soilstations.csv']);
listcell = textscan(fid, formatstr,'Headerlines', 1, 'Delimiter', ',');
fclose(fid);

sitesarray = [listcell{1}, listcell{2}];
clear listcell;

% Filter out sites for which there is no data in the data folder
havedatatest = ismember(sitesarray(:, 1), havedata);
sitesarray = sitesarray(havedatatest, :);

% Remove bad sites - only 583 so far
badsites = unique(badDataYears(:, 1));
badsitetest = ismember(sitesarray(:,1), badsites);
sitesarray(badsitetest, :) = [];

% Load data for each site into a large cellarray
datacell = cell(length(sitesarray), 1);
for i = 1:length(sitesarray);
    m = loadsnotel('daily', sitesarray(i,1));
    dailysite_datevec = datevec(m{2}, 'yyyy-mm-dd');
    % Create logical test for July and selected month
    julytest = dailysite_datevec(:,2)==7;
    monthtest = dailysite_datevec(:,2)==month;
    % Append mean July and selected month airtemps on to sitesarray
    sitesarray(i, 3) = nanmean(m{7}(julytest));% July mean max air temp
    sitesarray(i, 4) = nanmean(m{9}(julytest));% july mean avg air temp
    sitesarray(i, 5) = nanmean(m{15}(julytest));% july mean 20cm soil temp
    sitesarray(i, 6) = nanmean(m{7}(monthtest));% mean max air temp
    sitesarray(i, 7) = nanmean(m{9}(monthtest));% mean avg air temp
    sitesarray(i, 8) = nanmean(m{15}(monthtest));% mean 20cm soil temp
    sitesarray(i, 9) = nanmean(m{14}(monthtest));% mean 5cm soil temp
    % test if there is an average
    avgtest = averages(:,1)==sitesarray(i,1);
    if sum(avgtest)==1
        sitesarray(i, 10) = averages(avgtest,2); %AVG precip
        sitesarray(i, 11) = averages(avgtest,3); %AVG peakswe
    else
        sitesarray(i, 10) = NaN; %AVG precip
        sitesarray(i, 11) = nan; %AVG peakswe
    end
    sitesarray(i, 12) = (nanstd(m{9}(monthtest))./(nanmean(m{9}(monthtest))));% cv avg air temp
    sitesarray(i, 13) = nanstd(m{15}(monthtest));% cv 20cm soil temp
end

sitesarray(:,2) = (sitesarray(:,2)*.3048);


% PLOTS

% Full utah gradients
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name','Mean air and soil temps by elevation, Utah SNOTEL sites');

subplot (2, 2, 1)
plot(sitesarray(:,2), sitesarray(:,4), '.m');
hold on
plot(sitesarray(:,2), sitesarray(:,5), '.b');
xlabel('elevation(m)');
ylabel('Mean July temp (Celsius)');
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
title([ monthLabels{month} ' - Air and soil temp vs elevation']);
legend('Air', 'Soil(20cm)', 'Moist adiabatic lapse( 5^oC/km)');

subplot (2, 2, 3)
plot(sitesarray(:,2), sitesarray(:,3), '.r');
hold on
plot(sitesarray(:,2), sitesarray(:,6), '.b');
xlabel('elevation(m)');
ylabel('MAX air temp (Celsius)');
title('Max air temp vs elevation');
legend('July', [ monthLabels{month}, '(20cm, 12am)']);

subplot (2, 2, 4)
plot(sitesarray(:,2), sitesarray(:,10), '.k');
hold on
plot(sitesarray(:,2), sitesarray(:,11), '.b');
xlabel('elevation(m)');
ylabel('mm of water (Celsius)');
title('Precip and SWE');
legend('Total precip', 'Peak SWE');

% Restricted gradients
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name','Restricted gradients');

%Select peak SWE gradient
swetest = sitesarray(:, 11)<600 & sitesarray(:, 10)>400;
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


% Air temp gradient
fignum = fignum+1;    
h = figure(fignum);
set(h, 'Name', [ monthLabels{month} ' air temp vs elevation']);

plot(sitesarray(:,2), sitesarray(:,7), '.m');
hold on
% Plot moist adiabatic lapse rate
plot([1200 3500], [-1, -8.2], ':k')
% plot(sitesarray(:,2), sitesarray(:,7), '.b');
xlabel('elevation(m)');
ylabel('Mean air temp (Celsius)');
title([monthLabels{month} ' - Air temp vs elevation']);


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



fignum = fignum+1;    
h = figure(fignum);

badtest = sitesarray(:, 12)>5;

subplot (1,2,1);
errorbar(sitesarray(~badtest,2), sitesarray(~badtest,7), sitesarray(~badtest,12), 'ok', ...
    'MarkerFaceColor', [0.7 0.7 0.7]);
hold on
%plot(sitesarray(:,2), sitesarray(:,7), 'ok', 'MarkerFaceColor', ...
%    [0.7 0.7 0.7]);
errorbar(sitesarray(:,2), sitesarray(:,8), sitesarray(:,13), ...
    'ok', 'MarkerFaceColor', 'r');
%plot(sitesarray(:,2), sitesarray(:,8), 'ok', 'MarkerFaceColor', 'r');
% Plot moist adiabatic lapse rate
plot([900 3500], [4, -9], ':k')
xlabel('Elevation(m)');
ylabel('Mean temp (^oC)');
legend('Air', 'Soil (20cm)', 'Moist adiabatic lapse( 5^oC/km)');
title([ monthLabels{month} ' air & soil temp.']);

subplot(1,2,2);
offset = abs(sitesarray(:,7) - sitesarray(:,8));
plot(sitesarray(:,2), offset(:,1), 'ok', 'MarkerFaceColor', 'b');
xlabel('Elevation (m)');
ylabel('^oC');
title('Air-Soil temp. offset');
%legend('Air-Soil Offset');

junk = 99;
