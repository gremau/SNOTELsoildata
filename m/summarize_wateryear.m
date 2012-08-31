% summarize_wateryear.m
%
% Takes site filelists (from allsensors_daily, soilsensors_hourly),
% avgSWE, avgPrecip, and site inventory data and generates a series of
% files that contain climate or soil metrics that are aggregated for each
% wateryear present in the dataset (including a breakdown of monthly data
% and stats like elevation, avgSWE, avgPrecip, etc)
%
% Feb 20, 2012

clear;          % clear memory
close all;      % clear any figures
fignum=0;       % used to increment figure number for plots

addpath('/home/greg/data/programming_resources/m_common/'); % access to nanmean, etc

% Set data path and file name, read in file
rawdatapath = '../rawdata/';
processeddatapath = '../processed_data/';

% Other useful variables
wymonths = [10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9];

% Load list of sites with data in the daily data directory
allsiteslist = sortrows(dlmread([rawdatapath ...
    'allsensors_daily/sitelist.txt']));
soilsiteslist = sortrows(dlmread([rawdatapath ...
    'soilsensors_hourly/sitelist.txt']));

% Get 30yr average data
avgSWE = load7100Avg('swe');
avgPrecip = load7100Avg('precip');

% Get SNOTEL inventory data (elevation, coordinates)
siteinventory = sortrows(dlmread([processeddatapath 'SNOTELinventory.csv'], ...
    ',', 1, 5));

% Create three matrices to fill with climate and soil data
climatesummary = nan * zeros(length(allsiteslist), 42);
climatesummary(:, 1:2) = allsiteslist;

soilwatersummary = nan * zeros(length(soilsiteslist), 102);
soilwatersummary(:, 1:2) = soilsiteslist;

soiltempsummary = nan * zeros(length(soilsiteslist), 119);
soiltempsummary(:, 1:2) = soilsiteslist;


% CLIMATE SUMMARY (calculate, load, and write climatesummary)
site = 0;
for i = 1:length(climatesummary)
    % Load the daily datafile for the site (once per site - use for
    % multiple loops through the wateryears)
    if site~=climatesummary(i,1)
        site = climatesummary(i, 1);
        dailydata = loadsnotel('daily', site);
        sitedatevec = datevec(dailydata{2}, 'yyyy-mm-dd');
        % Get elevation, lat, lon, from siteinventory matrix
        inventoryrow = siteinventory(siteinventory(:, 1)==site, :);
        % Get Avg SWE data from avgSWE matrix
        avgSWErow = avgSWE(avgSWE(:, 1)==site, :);
        if isempty(avgSWErow)
            avgSWEvalue = nan;
        else
            avgSWEvalue = avgSWErow(27); % 30yr mean peak SWE
        end
        % Get Avg Precip data from avgPrecip matrix
        avgPreciprow = avgPrecip(avgPrecip(:, 1)==site, :);
        if isempty(avgPreciprow)
            avgPrecipvalue = nan;
        else
            avgPrecipvalue = avgPreciprow(15); % 30yr mean annual precip
        end
    end
    % Select the wateryear of interest
    wyear = climatesummary(i, 2);
    wyeartest = dailydata{21}==wyear; % use daily data wateryear
    wyeardatevec = sitedatevec(wyeartest, :);
    % Get peak swe and day of peak swe for the wateryear
    [maxswe, maxsweday] = nanmax(dailydata{4}(wyeartest));
    climatesummary(i, 3) = maxswe;
    if isnan(maxswe)
        maxsweday = nan; % make sure maxsweday is nan if maxswe also nan
    end
    climatesummary(i, 4) = maxsweday;
    climatesummary(i, 5) = nanmax(dailydata{10}(wyeartest)); % max snowdepth
    % Get all days that swe < 0.1, separate into early/late season days
    snowfree = dailydata{4}(wyeartest) < 0.1;
    snowfreedays = find(snowfree);
    earlyindex = snowfreedays < maxsweday;
    lateindex = snowfreedays > maxsweday;
    % Then find snowpack onset day, snowmelt day, snowpack duration,
    % airtemps during transitions
    if isempty(snowfreedays) % if no days are snow-free, set nans
        snowonsetday = nan;
        snowonsetdoy = nan;
        snowmeltday = nan;
        snowmeltdoy = nan; 
        snowpackduration = nan;
    elseif sum(earlyindex)==0 % no early-season snowfree days, onset = 1
        snowonsetdoy = 1;
        snowonsetday = datenum(wyeardatevec(1, :));
        snowmeltdoy = nanmin(snowfreedays(lateindex));
        snowmeltday = datenum(wyeardatevec(snowmeltdoy, :));
        snowpackduration = snowmeltday - snowonsetday;
    else
        snowonsetdoy = nanmax(snowfreedays(earlyindex));
        snowonsetday = datenum(wyeardatevec(snowonsetdoy, :));
        snowmeltdoy = nanmin(snowfreedays(lateindex));
        snowmeltday = datenum(wyeardatevec(snowmeltdoy, :));
        snowpackduration = snowmeltday - snowonsetday;
    end
    climatesummary(i, 6) = snowonsetday;
    climatesummary(i, 7) = snowmeltday;
    climatesummary(i, 8) = snowonsetdoy;
    climatesummary(i, 9) = snowmeltdoy;
    climatesummary(i, 10) = snowpackduration;
    % Accumulated and summer precip
    precip = dailydata{5}(wyeartest);
    accumprecip = precip(length(precip));
    climatesummary(i, 11) = accumprecip; % accumulated precip
    climatesummary(i, 12) = accumprecip - precip(182); % JAS precip
    % Monthly average air temps
    colindex1 = [13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35];
    % Get TAVG from data
    avgAirT = dailydata{9}(wyeartest);
    wydatevec = sitedatevec(wyeartest, :);
    for j = 1:12;
        monthtest = wydatevec(:, 2)==wymonths(j);
        meantemp = nanmean(avgAirT(monthtest));
        tempstd = nanstd(avgAirT(monthtest));
        climatesummary(i, colindex1(j)) = meantemp;
        climatesummary(i, colindex1(j) + 1) = tempstd;
    end
    climatesummary(i, 37) = nanmean(avgAirT); %MAT
    % Assign elev, lat, lon, avgSWE, avgPrecip from values set above
    climatesummary(i, 38) = inventoryrow(4)*0.3048; % elev
    climatesummary(i, 39) = inventoryrow(2); % lat
    climatesummary(i, 40) = inventoryrow(3); % lon
    climatesummary(i, 41) = avgSWEvalue;
    climatesummary(i, 42) = avgPrecipvalue;
end
    
clear dailydata wyear wyeartest monthtest;

% Write the climate data file
csvwrite([processeddatapath 'wyear_climatesummary.csv'], climatesummary);
    
site = 0;
for i = 1:length(soilsiteslist)
    % Load the daily AND hourly datafiles, plus datevecs and snowcover for
    % the site (use each in multiple loops through the wateryears)
    if site~=soilsiteslist(i, 1)
        site = soilsiteslist(i, 1);
        dailydata = loadsnotel('daily', site);
        hourlydata = loadsnotel('hourly', site);
        dailydatevec = datevec(dailydata{2}, 'yyyy-mm-dd');
        dailydatenum = datenum(dailydatevec);
        hourlydatevec = datevec(strcat(hourlydata{2}, hourlydata{3}), 'yyyy-mm-ddHH:MM');
        hourlydatenum = datenum(hourlydatevec);
        hourlysnowtest = swe_snowcover(site, hourlydatenum);
    end
    % Select the wateryear of interest
    wyear = soilsiteslist(i, 2);
    wyeartest_d = dailydata{21}==wyear;
    wyeartest_h = hourlydata{10}==wyear;
    %
    % SOIL WATER SUMMARY (calculate, load, and write soilwatersummary)
    %
    % Monthly soil moisture
    colindex2 = [3, 9, 15, 21, 27, 33, 39, 45, 51, 57, 61, 67];
    % Get soil moisture for 3 depths
%     sm5 = dailydata{11}(wyeartest_d);
%     sm20 = dailydata{12}(wyeartest_d);
%     sm50 = dailydata{13}(wyeartest_d);
    % NORMALIZE - Get soil moisture for 3 depths
    sm5 = smnormalize(dailydata{11}(wyeartest_d), 1);
    sm20 = smnormalize(dailydata{12}(wyeartest_d), 1);
    sm50 = smnormalize(dailydata{13}(wyeartest_d), 1);
    % Parse out the datevector for year i
    wydatevec_d = dailydatevec(wyeartest_d, :);
    % Loop through 12 months and average soil moisture for each
    for j = 1:12;
        monthtest_d = wydatevec_d(:, 2)==wymonths(j);
        mean5cm = nanmean(sm5(monthtest_d));
        std5cm = nanstd(sm5(monthtest_d));
        mean20cm = nanmean(sm20(monthtest_d));
        std20cm = nanstd(sm20(monthtest_d));
        mean50cm = nanmean(sm50(monthtest_d));
        std50cm = nanstd(sm50(monthtest_d));
        soilwatersummary(i, colindex2(j)) = mean5cm;
        soilwatersummary(i, colindex2(j) + 1) = std5cm;
        soilwatersummary(i, colindex2(j) + 2) = mean20cm;
        soilwatersummary(i, colindex2(j) + 3) = std20cm;
        soilwatersummary(i, colindex2(j) + 4) = mean50cm;
        soilwatersummary(i, colindex2(j) + 5) = std50cm;
    end
    % Create quarterly soil moisture tests
    ONDtest = (wydatevec_d(:, 2)==10) | (wydatevec_d(:, 2)==11) | ...
        (wydatevec_d(:, 2)==12);
    JFMtest = (wydatevec_d(:, 2)==1) | (wydatevec_d(:, 2)==2) | ...
        (wydatevec_d(:, 2)==3);
    AMJtest = (wydatevec_d(:, 2)==4) | (wydatevec_d(:, 2)==5) | ...
        (wydatevec_d(:, 2)==6);
    JAStest = (wydatevec_d(:, 2)==7) | (wydatevec_d(:, 2)==8) | ...
        (wydatevec_d(:, 2)==9);
    
    quarterscell = {ONDtest, JFMtest, AMJtest, JAStest};
    colindex3 = [73, 79, 85, 91];
    % Then calculate them and  put in soilwatersummary
    for j = 1:4;
        mean5cm = nanmean(sm5(quarterscell{j}));
        std5cm = nanstd(sm5(quarterscell{j}));
        mean20cm = nanmean(sm20(quarterscell{j}));
        std20cm = nanstd(sm20(quarterscell{j}));
        mean50cm = nanmean(sm50(quarterscell{j}));
        std50cm = nanstd(sm50(quarterscell{j}));
        soilwatersummary(i, colindex3(j)) = mean5cm;
        soilwatersummary(i, colindex3(j) + 1) = std5cm;
        soilwatersummary(i, colindex3(j) + 2) = mean20cm;
        soilwatersummary(i, colindex3(j) + 3) = std20cm;
        soilwatersummary(i, colindex3(j) + 4) = mean50cm;
        soilwatersummary(i, colindex3(j) + 5) = std50cm;
    end
    % Find the onset and snowmelt day using climatesummary
    findrow = (climatesummary(:, 1)==site & climatesummary(:, 2)==wyear);
    wyonsetday = climatesummary(findrow, 6);
    wymeltday = climatesummary(findrow, 7);
    % Calculate airT and soil moisture during the two-week transitions 
    % before onset, and before/after melt
    if isnan(wyonsetday)
        preonsetAirT = nan;
        preonset5cmSM = nan;
        preonset20cmSM = nan;
        preonset50cmSM = nan;
    else
        preonsettest = (dailydatenum<=wyonsetday) & ...
            (dailydatenum>=(wyonsetday-14)); % index the pre-onset interval
        preonsetAirT = nanmean(dailydata{9}(preonsettest));
        preonset5cmSM = nanmean(dailydata{11}(preonsettest));
        preonset20cmSM = nanmean(dailydata{12}(preonsettest));
        preonset50cmSM = nanmean(dailydata{13}(preonsettest));
    end
    % Same operations for AirT before and after melt
    if isnan(wymeltday)
        premeltAirT = nan;
        postmeltAirT = nan;
    else
        premelttest = (dailydatenum<=wymeltday & ...
            dailydatenum>=(wymeltday-14)); % index the premelt interval
        postmelttest = (dailydatenum>=wymeltday & ...
            dailydatenum<=(wymeltday+14)); % index the postmelt interval
        premeltAirT = nanmean(dailydata{9}(premelttest));
        postmeltAirT = nanmean(dailydata{9}(postmelttest));
    end
    % Put in the matrix
    soilwatersummary(i, 97) = preonsetAirT;
    soilwatersummary(i, 98) = preonset5cmSM;
    soilwatersummary(i, 99) = preonset20cmSM;
    soilwatersummary(i, 100) = preonset50cmSM;
    soilwatersummary(i, 101) = premeltAirT;
    soilwatersummary(i, 102) = postmeltAirT;
    %
    % SOIL TEMP SUMMARY (calculate, load, and write soiltempsummary
    %
    % Monthly soil temperature
    colindex4 = [3, 9, 15, 21, 27, 33, 39, 45, 51, 57, 61, 67];
    % Get soil temp for 3 depths
    st5 = hourlydata{7}(wyeartest_h);
    st20 = hourlydata{8}(wyeartest_h);
    st50 = hourlydata{9}(wyeartest_h);
    % Parse out the datevector and snowcover for year i
    wydatevec_h = hourlydatevec(wyeartest_h, :);
    wyhourlysnowtest = hourlysnowtest(wyeartest_h);
    % Loop through 12 months and average soil temp for each
    for j = 1:12;
        monthtest_h = wydatevec_h(:, 2)==wymonths(j);
        mean5cm = nanmean(st5(monthtest_h));
        std5cm = nanstd(st5(monthtest_h));
        mean20cm = nanmean(st20(monthtest_h));
        std20cm = nanstd(st20(monthtest_h));
        mean50cm = nanmean(st50(monthtest_h));
        std50cm = nanstd(st50(monthtest_h));
        soiltempsummary(i, colindex4(j)) = mean5cm;
        soiltempsummary(i, colindex4(j) + 1) = std5cm;
        soiltempsummary(i, colindex4(j) + 2) = mean20cm;
        soiltempsummary(i, colindex4(j) + 3) = std20cm;
        soiltempsummary(i, colindex4(j) + 4) = mean50cm;
        soiltempsummary(i, colindex4(j) + 5) = std50cm;
    end
    % Create quarterly soil temp tests
    ONDtest = (wydatevec_h(:, 2)==10) | (wydatevec_h(:, 2)==11) | ...
        (wydatevec_h(:, 2)==12);
    JFMtest = (wydatevec_h(:, 2)==1) | (wydatevec_h(:, 2)==2) | ...
        (wydatevec_h(:, 2)==3);
    AMJtest = (wydatevec_h(:, 2)==4) | (wydatevec_h(:, 2)==5) | ...
        (wydatevec_h(:, 2)==6);
    JAStest = (wydatevec_h(:, 2)==7) | (wydatevec_h(:, 2)==8) | ...
        (wydatevec_h(:, 2)==9);
    
    quarterscell = {ONDtest, JFMtest, AMJtest, JAStest};
    colindex5 = [73, 79, 85, 91];
    % Then calculate them and  put in soiltempsummary
    for j = 1:4;
        mean5cm = nanmean(st5(quarterscell{j}));
        std5cm = nanstd(st5(quarterscell{j}));
        mean20cm = nanmean(st20(quarterscell{j}));
        std20cm = nanstd(st20(quarterscell{j}));
        mean50cm = nanmean(st50(quarterscell{j}));
        std50cm = nanstd(st50(quarterscell{j}));
        soiltempsummary(i, colindex5(j)) = mean5cm;
        soiltempsummary(i, colindex5(j) + 1) = std5cm;
        soiltempsummary(i, colindex5(j) + 2) = mean20cm;
        soiltempsummary(i, colindex5(j) + 3) = std20cm;
        soiltempsummary(i, colindex5(j) + 4) = mean50cm;
        soiltempsummary(i, colindex5(j) + 5) = std50cm;
    end
    % Calculate some seasonal/yearly means
    soiltempsummary(i, 97) = nanmean(st5); % mast5cm
    soiltempsummary(i, 98) = nanstd(st5); % stdev
    soiltempsummary(i, 99) = nanmean(st20); % mast20cm
    soiltempsummary(i, 100) = nanstd(st20); % stdev
    soiltempsummary(i, 101) = nanmean(st50); % mast50cm
    soiltempsummary(i, 102) = nanstd(st50); % stdev
    % Snowcovered soil temp means
    soiltempsummary(i, 103) = nanmean(st5(wyhourlysnowtest));
    soiltempsummary(i, 104) = nanstd(st5(wyhourlysnowtest));
    soiltempsummary(i, 105) = nanmean(st20(wyhourlysnowtest));
    soiltempsummary(i, 106) = nanstd(st20(wyhourlysnowtest));
    soiltempsummary(i, 107) = nanmean(st50(wyhourlysnowtest));
    soiltempsummary(i, 108) = nanstd(st50(wyhourlysnowtest));
    % Snowfree soil temp means
    soiltempsummary(i, 109) = nanmean(st5(~wyhourlysnowtest));
    soiltempsummary(i, 110) = nanstd(st5(~wyhourlysnowtest));
    soiltempsummary(i, 111) = nanmean(st20(~wyhourlysnowtest));
    soiltempsummary(i, 112) = nanstd(st20(~wyhourlysnowtest));
    soiltempsummary(i, 113) = nanmean(st50(~wyhourlysnowtest));
    soiltempsummary(i, 114) = nanstd(st50(~wyhourlysnowtest));
    % Get soil temperatures at transition times (see above)
    wydays_h = find(wyeartest_h); % get full dataset indices for water year
    if isnan(wyonsetday)
        preonset5cmST = nan;
        preonset20cmST = nan;
        preonset50cmST = nan;
    else
        preonsettest = (hourlydatenum<=wyonsetday & ...
            hourlydatenum>=(wyonsetday-14)); % index the pre-onset interval
        preonset5cmST = nanmean(hourlydata{7}(preonsettest));
        preonset20cmST = nanmean(hourlydata{8}(preonsettest));
        preonset50cmST = nanmean(hourlydata{9}(preonsettest));
    end
    % Same operations for SoilT before and after melt
    if isnan(wymeltday)
        premelt5cmST = nan;
        postmelt5cmST = nan;
    else 
        premelttest = (hourlydatenum<=wymeltday & ...
            hourlydatenum>=(wymeltday-14)); % index the premelt interval
        postmelttest = (hourlydatenum>=wymeltday & ...
            hourlydatenum<=(wymeltday+14)); % index the postmelt interval
        premelt5cmST = nanmean(hourlydata{7}(premelttest));
        postmelt5cmST = nanmean(hourlydata{7}(postmelttest));
    end
    % Put in the matrix
    soiltempsummary(i, 115) = preonset5cmST;
    soiltempsummary(i, 116) = preonset20cmST;
    soiltempsummary(i, 117) = preonset50cmST;
    soiltempsummary(i, 118) = premelt5cmST;
    soiltempsummary(i, 119) = postmelt5cmST;
    
end

% Write the soil data files
% csvwrite([processeddatapath 'wyear_soilwatersummary.csv'], soilwatersummary);
csvwrite([processeddatapath 'wyear_soilwatersummary_smnorm.csv'], soilwatersummary);
csvwrite([processeddatapath 'wyear_soiltempsummary.csv'], soiltempsummary);

junk = 99;
