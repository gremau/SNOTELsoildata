% plot_monthly_ts_vs_elev.m
%
% ver 3: 111130 GM
%
%
% Was monthlysoiltempgradients_scr.m

close all;      % clear any figures
fignum = 0;     % used to increment figure number for plots

% Set data path and file name, read in file
datapath = '../rawdata/';
%datapath = '/home/greg/data/rawdata/SNOTELdata/';

% Generate list of sites from the _sitelist.txt file
haveData = unique(dlmread([datapath 'soilsensors_hourly/sitelist.txt']));

% Generate list of sites and their elevations from inventory file
% Create format string (station,elev only here)
formatstr = '%*s%f%*s%*s%*s%*s%f%*s%*s%*s%*s%*s%*s%*s%*s%*u';
fid = fopen([datapath 'station_inventory/UT_soilstations.csv']);
%fid = fopen([datapath 'station_invventory/merged_soilstations.csv']);
listcell = textscan(fid, formatstr,'Headerlines', 1, 'Delimiter', ',');
fclose(fid);

% Get 30yr average data
avgSWE = load7100Avg('swe');
avgPrecip = load7100Avg('precip');


sitesArray = [listcell{1}, listcell{2}];
precipArray = [precipcell{1}, precipcell{2}, precipcell{3}];
clear listcell;

% Filter out sites for which there is no data in the data folder
testHaveData = ismember(sitesArray(:, 1), haveData);
sitesArray = sitesArray(testHaveData, :);
testHaveData2 = ismember(precipArray(:, 1), haveData);
precipArray = precipArray(testHaveData2, :);

% Load list of bad data years for all sites
badDataYears = sortrows(csvread([datapath 'allsensors_daily/baddata.txt'], 1, 0));

% Change variables to something more useful
sites = sitesArray(:,1);
elev = sitesArray(:,2);
ltprecip = nan * zeros(length(sites), 1);
ltswe = nan * zeros(length(sites), 1);

for i = 1:length(sites);
    test = ismember(precipArray(:,1), sites(i));
    if sum(test)==1
        ltprecip(i) = precipArray(test, 2);
        ltswe(i) = precipArray(test, 3);
    else
        ltprecip(i) = nan;
        ltswe(i) = nan;
    end
end

t = zeros(13, length(sites));
template = {t t t t t t};
[snowmeans freemeans snowsums freesums snowdev freedev] = template{:};

m = cell(1, length(sites));

for i = 1:length(sites);
    [m, ~] = loadsnotel('hourly', sites(i));
    Ts = m{7};
    
    % INTERPOLATION to fill gaps in Ts (helps with running mean and
    % variance calculations below.
    %
    % Interpolates over NaNs (data gaps) in the input time series (may be
    % complex), but ignores trailing and leading NaN.
    %
    % from FIXGAPS routine on Matlab Central file exchange,
    % by R. Pawlowicz 6/Nov/99

    Ts_filled = Ts;

    bad = isnan(Ts);
    good = find(~bad);

    bad([1:(min(good)-1) (max(good)+1):end]) = 0;

    Ts_filled(bad)=interp1(good, Ts(good), find(bad), 'pchip');
    
    Ts = Ts_filled;
    %
    % Filter
    %Ts = filterseries(Ts, 'shift', 2.5);
    % FILTER Tsoil data (returns filtered and re-interpolated array) 
    Ts_meandiff = filterseries(Ts, 'mean', 7);
    Ts = filterseries(Ts_meandiff, 'shift', 5);
    
    % Generate decimal day and snowcover test arrays
    decday_h = datenum(strcat(m{2}, m{3}), 'yyyy-mm-ddHH:MM');
    date_vec = datevec(strcat(m{2}, m{3}), 'yyyy-mm-ddHH:MM');
    test_snow = swe_snowcover(sites(i), decday_h);
    test_free = ~test_snow;
    %nan_Ts = isnan(Ts); % Could be changed to nanmean
    
    % Plot filtered and filled Ts data and snowcover test (QC) 
%     h = figure;
%     set(h, 'Name', ['Site ' num2str(site_id) ' - Snow/snowfree mean Ts']);
%     subplot 311;
%     plot(decday_h, Ts, 'r');
%     subplot 312;
%     plot(decday_h, test_snow);
%     subplot 313;
%     plot(decday_h, test_free);
    
    for j = 1:12
        test_month = date_vec(:,2) == j;
        snowmeans(j, i) = mean(Ts(test_month & test_snow,:));
        snowdev(j, i) = std(Ts(test_month & test_snow,:));
        snowsums(j, i) = sum(Ts(test_month & test_snow,:));
        freemeans(j, i) = mean(Ts(test_month & test_free,:));
        freedev(j, i) = std(Ts(test_month & test_free,:));
        freesums(j, i) = sum(Ts(test_month & test_free,:));
    end
    
    snowmeans(13, i) = mean(Ts(test_snow,:));
    snowdev(13, i) = std(Ts(test_snow,:));
    snowsums(13, i) = std(Ts(test_snow,:));
    freemeans(13, i) = mean(Ts(test_free,:));
    freedev(13, i) = std(Ts(test_free,:));
    freesums(13, i) = sum(Ts(test_free,:));
    
    clear m decday_h test_snow test_free test_month
end

clear headers Ts nan_Ts i j;

% convert to m
elev = elev/3.28;

% % Plot all temperature data
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Mean 5cm Ts, by elevation - All months/years');

errorbar(elev, snowmeans(13,:), snowdev(13,:), 'b.');
hold on
errorbar(elev, freemeans(13,:), freedev(13,:), 'r.');
plot(elev, snowmeans(13,:), 'b.', elev, freemeans(13,:), 'r.')
title('Mean soil temperature by elevation');
legend('Snow-covered', 'Snow-free');
ylabel('Mean 5cm soil temp');
xlabel('Elevation (m)');

% Plot February and July soil temps by elevation
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Mean 5cm Ts, Feb and July, by elevation');

%subplot 221
errorbar(elev, snowmeans(2,:), snowdev(2,:), 'ko', 'MarkerFaceColor', 'w');
hold on;
errorbar(elev, freemeans(7,:), freedev(7,:), 'ko', 'MarkerFaceColor', 'k');
% Plot moist adiabatic lapse rate
plot([1600 3400], [18, 9], ':k');
legend('February', 'July', 'Moist adiabatic lapse (5^oC/km)');
xlabel('Elevation (m)');
ylabel('Mean soil temp (^oC)')


% Plot February and July soil temps by SWE
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Mean 5cm Ts, Feb and July, by mean MAX SWE');
errorbar(ltswe, snowmeans(2,:), snowdev(2,:), 'b.');
hold on
errorbar(ltswe, freemeans(7,:), freedev(7,:), 'r.');
plot(ltswe, snowmeans(2,:), 'b.', elev, freemeans(7,:), 'r.')
title('5cm Ts: Feb = blue, July = red');
xlabel('SWE (mm)');
ylabel('Mean degrees C')


% Plot by month
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Mean 5cm Ts, by elevation, by month');

for i = 1:12;
    subplot (3, 4, i);
    errorbar(elev, snowmeans(i,:), snowdev(i,:), 'b.');
    hold on
    errorbar(elev, freemeans(i,:), freedev(i,:), 'r.');
    plot(elev, snowmeans(i,:), 'b.', elev, freemeans(i,:), 'r.')
    title(num2str(i));
end

% Plot by month - Fixed axes
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Mean 5cm Ts, by elevation, by month (fixed axes)');

for i = 1:12;
    subplot (3, 4, i);
    errorbar(elev, snowmeans(i,:), snowdev(i,:), 'b.');
    hold on
    errorbar(elev, freemeans(i,:), freedev(i,:), 'r.');
    plot(elev, snowmeans(i,:), 'b.', elev, freemeans(i,:), 'r.')
    title(num2str(i));
    ylim([-5 30]);
end

% Plot summed Ts by month
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Summed 5cm Ts, by elevation, by month');

for i = 1:12;
    subplot (3, 4, i);
    plot(elev, snowsums(i,:), 'b.', elev, freesums(i,:), 'r.')
    title(num2str(i));
end  

junk = 99;
