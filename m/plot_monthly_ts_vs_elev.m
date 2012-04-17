% plot_monthly_ts_vs_elev.m
%
% ver 3: 111130 GM
%
%
% Was monthlysoiltempgradients_scr.m

close all;      % clear any figures
fignum = 0;     % used to increment figure number for plots

% Set data path and file name, read in file
rawdatapath = '../rawdata/';
processeddatapath = '../processed_data/';

% Ask user for month number and state
% monthsel = str2double(input('Which month (1-12)?: ', 's'));
statesel = input(...
    'Which state ("AZ, CO, ID, MT, NM, NV, UT, WY, or all")?: ', 's');

% monthlabels = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sept' 'Oct'...
%     'Nov' 'Dec'};
% monthlabel = monthlabels{monthsel};

%Load list of sites in the daily data directory
havedata = csvread([rawdatapath 'soilsensors_hourly/sitelist.txt']);
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

% Create some useful vectors
states = inventorycell{1};
sitesarray = [inventorycell{2}, inventorycell{5}*0.3048];
precipArray = [avgPrecip(:, 1), avgPrecip(:, 15)];
SWEarray = [avgSWE(:, 1), avgSWE(:, 27)];
clear inventorycell;

% Remove site rows that are not part of the selected state
if strcmpi(statesel, 'all')
    sitesarray = sitesarray;
else
    stateseltest = strcmpi(states, statesel);
    sitesarray = sitesarray(stateseltest, :);
end

% Filter out sites for which there is no data in the data folder
testhavedata = ismember(sitesarray(:, 1), havedata);
sitesarray = sitesarray(testhavedata, :);

% Load list of bad data years for all sites
%badDataYears = sortrows(csvread([rawdatapath 'allsensors_daily/baddata.txt'], 1, 0));

% Change variables to something more useful
% sites = sitesArray(:,1);
% elev = sitesArray(:,2);
% ltprecip = nan * zeros(length(sites), 1);
% ltswe = nan * zeros(length(sites), 1);

% Add 30-year SWE and Precip for all sites remaining in sitesarray
for i = 1:length(sitesarray(:,1));
    avgPreciptest = ismember(avgPrecip(:, 1), sitesarray(i, 1));
    avgSWEtest = ismember(avgSWE(:,1), sitesarray(i, 1));

    if sum(avgSWEtest)==1
        sitesarray(i, 3) = avgSWE(avgSWEtest, 27); % AVG peakswe
    else
        sitesarray(i, 3) = nan;
    end
    if sum(avgPreciptest)==1
        sitesarray(i, 4) = avgPrecip(avgPreciptest, 15); % AVG precip
    else
        sitesarray(i, 4) = nan;
    end
end

% Create 6 empty arrays to accumulate calculated values
t = zeros(13, length(sitesarray));
template = {t t t t t t};
[snowmeans freemeans snowsums freesums snowdev freedev] = template{:};

%m = cell(1, length(sitesarray));

for i = 1:length(sitesarray(:, 1));
    m = loadsnotel('hourly', sitesarray(i, 1));
    Ts = m{8}; % Select sensor depth (7=5cm, 8=20cm, 9=50cm)
    
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
    Ts_meandiff = filtertempseries(Ts, 'mean', 7);
    Ts = filtertempseries(Ts_meandiff, 'shift', 5);
    
    % Generate decimal day and snowcover test arrays
    decday_h = datenum(strcat(m{2}, m{3}), 'yyyy-mm-ddHH:MM');
    date_vec = datevec(strcat(m{2}, m{3}), 'yyyy-mm-ddHH:MM');
    test_snow = swe_snowcover(sitesarray(i, 1), decday_h);
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


% % Plot all temperature data
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', 'Mean 5cm Ts, by elevation - All months/years');

errorbar(sitesarray(:, 2) , snowmeans(13,:), snowdev(13,:), 'b.');
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
