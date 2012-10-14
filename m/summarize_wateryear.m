% summarize_wateryear.m
%
% Takes site filelists (from allsensors_daily, soilsensors_hourly),
% avgSWE, avgPrecip, and site inventory data and generates a series of
% files that contain climate or soil metrics that are aggregated for each
% wateryear present in the dataset (including a breakdown of monthly data
% and stats like elevation, avgSWE, avgPrecip, etc). Bad data are removed 
% and remaining data are filtered with a 3 sigma filter before calculation.
%
% Options: 
% 1: Use daily or hourly soil data (Ts and VWC)
% 2: Normalize soil moisture data (after filtering)
%
% Feb 20, 2012
%
clear;          % clear memory
close all;      % clear any figures
fignum=0;       % used to increment figure number for plots

% Ask user input
soilinput = input('Use hourly or daily soil data?: ', 's');
if strcmpi(soilinput, 'daily')
    windowsize = 11;
    filtlist = [4, 6:16]; % don't filter precip
    % Mark where the needed soil columns are in the datafile
    smcol = [11, 12, 13];
    stcol = [14, 15, 16];
elseif strcmpi(soilinput, 'hourly')
    windowsize = 25;
    filtlist = 4:9;
    smcol = [4, 5, 6];
    stcol = [7, 8, 9];
end
norminput = input('Normalize soil moisture data? (y/n): ', 's');
if strcmpi(norminput, 'y')
    normstr = '_smnorm';
else
    normstr = '';
end

% Get access to nanmean, etc
addpath('/home/greg/data/code_resources/m_common/nanstuff/'); 

% Set data path and file name, read in file
rawdatapath = '../rawdata/';
processeddatapath = '../processed_data/';

% Other useful variables
wymonths = [10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9];

% Load list of sites with data in the daily data directory
allsiteslist = sortrows(dlmread([rawdatapath ...
    'allsensors_daily/filelist.txt']));
soilsiteslist = sortrows(dlmread([rawdatapath ...
    'soilsensors_hourly/filelist.txt']));

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

% *** BEGIN CLIMATE CALCULATIONS ****************************************
% climatesummary matrix is generated in this loop, with one interation per 
% site/wateryear.

site = 0;
for i = 1:length(climatesummary)
    % Load the daily datafile for the site (once per site - use for
    % multiple loops through the wateryears)
    if site~=climatesummary(i,1)
        site = climatesummary(i, 1);
        dailydata = loadsnotel(site, 'daily', 'exclude');
        % Filter SWE, temp, and depth data. Windowsize=11 days
        for j = [4, 6:10]
            dailydata{j} = filterseries(dailydata{j}, 'sigma', 11, 3);
        end
        % Get a date vector and a datenum array
        datevec_c = datevec(dailydata{2}, 'yyyy-mm-dd');
        datenum_c = datenum(datevec_c);
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
    % Select the wateryear of interest and create a wateryear datvevec
    wyear_c = climatesummary(i, 2);
    wyeartest_c = dailydata{21}==wyear_c; % use daily data wateryear
    wyeardatevec_c = datevec_c(wyeartest_c, :);
    
    % SNOWPACK - Peak, duration, onset, etc
    %
    % Get peak swe and day of peak swe for the wateryear
    [maxswe, maxsweday] = nanmax(nancheck(dailydata{4}(wyeartest_c)));
    climatesummary(i, 3) = maxswe;
    if isnan(maxswe)
        maxsweday = nan; % make sure maxsweday is nan if maxswe also nan
    end
    climatesummary(i, 4) = maxsweday;
    % Max snow depth
    climatesummary(i, 5) = nanmax(nancheck(dailydata{10}(wyeartest_c)));
    % Then find snowpack onset day, snowmelt day, snowpack duration...
    % First, mark days where swe < 0.1, WARNING - NaN's are counted as snow
    snowfree = nancheck(dailydata{4}(wyeartest_c)) < 0.1;
    if sum(snowfree)==0 % if no days are snow-free, set nans
        totaldaysSC = nan; % Total SC days
        maxcontinSC = nan; % Length of longest snowcovered period
        numcontinSC = nan; % Number of snowcovered periods
        snowonsetday = nan;
        snowmeltday = nan;
        snowpackduration = nan;
    elseif sum(snowfree)>363 % if there are no snow days, set zeros, nans
        totaldaysSC = 0;
        maxcontinSC = 0;
        numcontinSC = 0;
        snowonsetday = nan;
        snowmeltday = nan;
        snowpackduration = 0;
    else
        snowfree(330:end) = 1; % Ignore late season snow events
        % Cumsum the boolean array
        cumSnowfree = cumsum(snowfree);
        % Put this in a frequency histogram to identify contiguous snowcovers
        xedges = 0:max(cumSnowfree)+1;
        snowhist = histc(cumSnowfree, xedges);
        snowperiods = snowhist(snowhist>2); % Ignore snow present < 1 day
        % Get some statistics
        totaldaysSC = length(cumSnowfree) - max(cumSnowfree); % Total SC days
        maxcontinSC = max(snowperiods); % Length of longest snowcovered period
        numcontinSC = length(snowperiods); % Number of snowcovered periods
        %Find the onset and snowmelt day using snowhist, cumSnowfree
        if snowfree(1)==0
            snowonsetday = min(find(cumSnowfree==...
                min(find(snowhist==snowperiods(1))-1)));
        else
            snowonsetday = min(find(cumSnowfree==...
                min(find(snowhist==snowperiods(1))-1)))+1;
        end
        snowmeltday = max(find(cumSnowfree==...
            max(find(snowhist==snowperiods(end)))));
    end
    climatesummary(i, 6) = snowonsetday;
    climatesummary(i, 7) = snowmeltday;
    climatesummary(i, 8) = snowmeltday-snowonsetday; % Duration of snowpack
    climatesummary(i, 9) = totaldaysSC; % Total days, may vary from above
    climatesummary(i, 10) = maxcontinSC;
    climatesummary(i, 11) = numcontinSC;
    
    % PRECIP - Accumulated and summer precip
    precip = dailydata{5}(wyeartest_c);
    accumprecip = precip(length(precip));
    climatesummary(i, 12) = accumprecip; % accumulated precip
    if length(precip) > 355 % Allow for 10 days missing precip data
        climatesummary(i, 13) = accumprecip - precip(182); % JAS precip
    else
        climatesummary(i, 13) = NaN;
    end
    
    % Monthly average SWE
    colindex1 = 14:3:47;
    % Get average SWE from data
    wteq = dailydata{4}(wyeartest_c);
    for j = 1:12;
        monthtest = wyeardatevec_c(:, 2)==wymonths(j);
        meanswe = nanmean(nancheck(wteq(monthtest)));
        slice = nancheck(wteq(monthtest));
        medswe = median(slice(~isnan(slice)));
        swestd = nanstd(nancheck(wteq(monthtest)));
        climatesummary(i, colindex1(j)) = meanswe;
        climatesummary(i, colindex1(j) + 1) = medswe;
        climatesummary(i, colindex1(j) + 2) = swestd;
    end
    clear monthtest;
    
    % Monthly average air temps
    colindex2 = 50:2:72;
    % Get TAVG from data
    avgAirT = dailydata{9}(wyeartest_c);
    for j = 1:12;
        monthtest = wyeardatevec_c(:, 2)==wymonths(j);
        meantemp = nanmean(nancheck(avgAirT(monthtest)));
        tempstd = nanstd(nancheck(avgAirT(monthtest)));
        climatesummary(i, colindex2(j)) = meantemp;
        climatesummary(i, colindex2(j) + 1) = tempstd;
    end
    clear monthtest;
    climatesummary(i, 74) = nanmean(nancheck(avgAirT)); % MAT
    climatesummary(i, 75) = nanstd(nancheck(avgAirT)); % MAT - sd
    
    % Calculate air temp during the two-week transitions 
    % before onset of snowpack an at snowmelt
    %
    % Get the datenums of melt day and onset day  
    minwydatenum_c = min(floor(datenum_c(wyeartest_c)));
    onsetdatenum_c = snowonsetday + minwydatenum_c;
    meltdatenum_c = snowmeltday + minwydatenum_c;
    % Find the mean airT for the 2 weeks prior to snow onset
    if isnan(onsetdatenum_c)
        preonsetAirT = nan;
        preonsetAirTsd = nan;
    else
        preonsettest = (datenum_c<=onsetdatenum_c) & ...
            (datenum_c>=(onsetdatenum_c-14)); % index the pre-onset interval
        preonsetAirT = nanmean(nancheck(dailydata{9}(preonsettest)));
        preonsetAirTsd = nanstd(nancheck(dailydata{9}(preonsettest)));
    end
    
    % Same operations for AirT before and after melt
    if isnan(meltdatenum_c)
        premeltAirT = nan;
        premeltAirTsd = nan;
        postmeltAirT = nan;
        postmeltAirTsd = nan;
    else
        premelttest = (datenum_c<=meltdatenum_c & ...
            datenum_c>=(meltdatenum_c-14)); % index the premelt interval
        postmelttest = (datenum_c>=meltdatenum_c & ...
            datenum_c<=(meltdatenum_c+14)); % index the postmelt interval
        premeltAirT = nanmean(nancheck(dailydata{9}(premelttest)));
        premeltAirTsd = nanstd(nancheck(dailydata{9}(premelttest)));
        postmeltAirT = nanmean(nancheck(dailydata{9}(postmelttest)));
        postmeltAirTsd = nanstd(nancheck(dailydata{9}(postmelttest)));
    end
    climatesummary(i, 76) = preonsetAirT;
    climatesummary(i, 77) = preonsetAirTsd;
    climatesummary(i, 78) = premeltAirT;
    climatesummary(i, 79) = premeltAirTsd;
    climatesummary(i, 80) = postmeltAirT;
    climatesummary(i, 81) = postmeltAirTsd;
    % Assign elev, lat, lon, avgSWE, avgPrecip from values set above
    climatesummary(i, 82) = inventoryrow(4)*0.3048; % elev
    climatesummary(i, 83) = inventoryrow(2); % lat
    climatesummary(i, 84) = inventoryrow(3); % lon
    climatesummary(i, 85) = avgSWEvalue;
    climatesummary(i, 86) = avgPrecipvalue;
end
    
clear i j dailydata wyear_c wyeartest_c monthtest slice;

% Write the climate data file
csvwrite([processeddatapath 'wyear_climatesummary.txt'], climatesummary);


% *** BEGIN SOIL CALCULATIONS ******************************************
% Both soiltempsummary and soilwatersummary are generated in one loop, with
% one interation per site/wateryear.
% IF soilinput is 'daily', daily data is loaded, hourly if 'hourly'
% IF norminput is y, soil moisture data are normalized

site = 0;
for i = 1:length(soilsiteslist)
    % When we loop into a new site in soilsitelist:
    % Load the datafile, plus datevecs and snowcover for the site (use 
    % each in multiple loops through the wateryears)
    if site~=soilsiteslist(i, 1)
        site = soilsiteslist(i, 1);
        % Load daily or hourly data based on user input
        data = loadsnotel(site, soilinput, 'exclude');
        % Filter data using filtlist
        % Daily = SWE, temp, depth, and VWC data. windowsize=11 days
        % Hourly = all soil sensors, windowsize = 25
        for j = filtlist
            data{j} = filterseries(data{j}, 'sigma', windowsize, 3);
        end
        % Normalize, based on user input
        if strcmpi(norminput, 'y')
            for j = smcol
                data{j} = smnormalize(data{j}, 1);
            end
        end
        % Get datevecs and datenum arrays, depending on user input
        if strcmpi(soilinput, 'daily')
            datevec_s = datevec(data{2}, 'yyyy-mm-dd');
            datenum_s = datenum(datevec_s);
            snowtest = data{4} > 0.1;
            wyearvec_s = data{21};
        elseif strcmpi(soilinput, 'hourly')
            datevec_s = datevec(strcat(data{2}, data{3}), 'yyyy-mm-ddHH:MM');
            datenum_s = datenum(datevec_s);
            snowtest = swe_snowcover(site, datenum_s);
            wyearvec_s = data{10};
        end
    end
    % Select wateryear i and create a datevector for wateryear i
    wyear_s = soilsiteslist(i, 2);
    wyeartest_s = wyearvec_s==wyear_s;
    wydatevec_s = datevec_s(wyeartest_s, :);
    
    % SOIL WATER SUMMARY (calculate, load, and write soilwatersummary)
    %
    % Monthly soil moisture column index
    colindex3 = [3, 9, 15, 21, 27, 33, 39, 45, 51, 57, 63, 69];
    % Assign the data for wateryear i
    % NOTE - normalizing each year individually
    % is also an option
    sm5  = data{smcol(1)}(wyeartest_s);
    sm20 = data{smcol(2)}(wyeartest_s);
    sm50 = data{smcol(3)}(wyeartest_s);
%     if strcmpi(norminput, 'n')
%         % Assign soil moisture data for 3 depths
%         sm5  = data{smcol(1)}(wyeartest_s);
%         sm20 = data{smcol(2)}(wyeartest_s);
%         sm50 = data{smcol(3)}(wyeartest_s);
%     elseif strcmpi(norminput, 'y')
%         % normalize with smnorm
%         sm5  = smnormalize(data{smcol(1)}(wyeartest_s), 1);
%         sm20 = smnormalize(data{smcol(2)}(wyeartest_s), 1);
%         sm50 = smnormalize(data{smcol(3)}(wyeartest_s), 1);
%     end

    % Loop through 12 months and average soil moisture for each
    for j = 1:12;
        monthtest = wydatevec_s(:, 2)==wymonths(j);
        mean5cm = nanmean(nancheck(sm5(monthtest)));
        std5cm = nanstd(nancheck(sm5(monthtest)));
        mean20cm = nanmean(nancheck(sm20(monthtest)));
        std20cm = nanstd(nancheck(sm20(monthtest)));
        mean50cm = nanmean(nancheck(sm50(monthtest)));
        std50cm = nanstd(nancheck(sm50(monthtest)));
        soilwatersummary(i, colindex3(j)) = mean5cm;
        soilwatersummary(i, colindex3(j) + 1) = std5cm;
        soilwatersummary(i, colindex3(j) + 2) = mean20cm;
        soilwatersummary(i, colindex3(j) + 3) = std20cm;
        soilwatersummary(i, colindex3(j) + 4) = mean50cm;
        soilwatersummary(i, colindex3(j) + 5) = std50cm;
    end
    clear monthtest;
    
    % Create quarterly soil moisture tests
    ONDtest = (wydatevec_s(:, 2)==10) | (wydatevec_s(:, 2)==11) | ...
        (wydatevec_s(:, 2)==12);
    JFMtest = (wydatevec_s(:, 2)==1) | (wydatevec_s(:, 2)==2) | ...
        (wydatevec_s(:, 2)==3);
    AMJtest = (wydatevec_s(:, 2)==4) | (wydatevec_s(:, 2)==5) | ...
        (wydatevec_s(:, 2)==6);
    JAStest = (wydatevec_s(:, 2)==7) | (wydatevec_s(:, 2)==8) | ...
        (wydatevec_s(:, 2)==9);
    
    quarterscell = {ONDtest, JFMtest, AMJtest, JAStest};
    colindex4 = [75, 81, 87, 93];
    % Then calculate them and  put in soilwatersummary
    for j = 1:4;
        mean5cm = nanmean(nancheck(sm5(quarterscell{j})));
        std5cm = nanstd(nancheck(sm5(quarterscell{j})));
        mean20cm = nanmean(nancheck(sm20(quarterscell{j})));
        std20cm = nanstd(nancheck(sm20(quarterscell{j})));
        mean50cm = nanmean(nancheck(sm50(quarterscell{j})));
        std50cm = nanstd(nancheck(sm50(quarterscell{j})));
        soilwatersummary(i, colindex4(j)) = mean5cm;
        soilwatersummary(i, colindex4(j) + 1) = std5cm;
        soilwatersummary(i, colindex4(j) + 2) = mean20cm;
        soilwatersummary(i, colindex4(j) + 3) = std20cm;
        soilwatersummary(i, colindex4(j) + 4) = mean50cm;
        soilwatersummary(i, colindex4(j) + 5) = std50cm;
    end
    
    % Calculate soil moisture during the two-weeks prior to snowpack onset
    %
    % First get the datenums of melt day and onset day
    findrow = (climatesummary(:, 1)==site & climatesummary(:, 2)==wyear_s);   
    minwydatenum_s = min(floor(datenum_s(wyeartest_s)));
    onsetdatenum_s = climatesummary(findrow, 6) + minwydatenum_s;
    meltdatenum_s = climatesummary(findrow, 7) + minwydatenum_s;
    % Calculate soil moisture during the two-week transitions 
    % before onset, and before/after melt
    if isnan(onsetdatenum_s)
        preonsetAirT = nan;
        preonsetAirTsd = nan;
        preonset5cmSM = nan;
        preonset5cmSMsd = nan;
        preonset20cmSM = nan;
        preonset20cmSMsd = nan;
        preonset50cmSM = nan;
        preonset50cmSMsd = nan;
    else
        preonsettest = (datenum_s<=onsetdatenum_s) & ...
            (datenum_s>=(onsetdatenum_s-14)); % index the pre-onset interval
        preonset5cmSM = nanmean(nancheck(data{smcol(1)}(preonsettest)));
        preonset5cmSMsd = nanstd(nancheck(data{smcol(1)}(preonsettest)));
        preonset20cmSM = nanmean(nancheck(data{smcol(2)}(preonsettest)));
        preonset20cmSMsd = nanstd(nancheck(data{smcol(2)}(preonsettest)));
        preonset50cmSM = nanmean(nancheck(data{smcol(3)}(preonsettest)));
        preonset50cmSMsd = nanstd(nancheck(data{smcol(3)}(preonsettest)));
    end

    % Put in the matrix
    soilwatersummary(i, 99) = preonset5cmSM;
    soilwatersummary(i, 100) = preonset5cmSMsd;
    soilwatersummary(i, 101) = preonset20cmSM;
    soilwatersummary(i, 102) = preonset20cmSMsd;
    soilwatersummary(i, 103) = preonset50cmSM;
    soilwatersummary(i, 104) = preonset50cmSMsd;

    % SOIL TEMP SUMMARY (calculate, load, and write soiltempsummary
    %
    % Monthly soil temperature column indices
    colindex5 = [3, 9, 15, 21, 27, 33, 39, 45, 51, 57, 63, 69];
    % Get wateryear soil temp for 3 depths
    st5  = data{stcol(1)}(wyeartest_s);
    st20 = data{stcol(2)}(wyeartest_s);
    st50 = data{stcol(3)}(wyeartest_s);
    
    % Parse out the snowcover for year i
    wysnowtest = snowtest(wyeartest_s);
    
    % Loop through 12 months and average soil temp for each
    for j = 1:12;
        monthtest = wydatevec_s(:, 2)==wymonths(j);
        mean5cm = nanmean(nancheck(st5(monthtest)));
        std5cm = nanstd(nancheck(st5(monthtest)));
        mean20cm = nanmean(nancheck(st20(monthtest)));
        std20cm = nanstd(nancheck(st20(monthtest)));
        mean50cm = nanmean(nancheck(st50(monthtest)));
        std50cm = nanstd(nancheck(st50(monthtest)));
        soiltempsummary(i, colindex5(j)) = mean5cm;
        soiltempsummary(i, colindex5(j) + 1) = std5cm;
        soiltempsummary(i, colindex5(j) + 2) = mean20cm;
        soiltempsummary(i, colindex5(j) + 3) = std20cm;
        soiltempsummary(i, colindex5(j) + 4) = mean50cm;
        soiltempsummary(i, colindex5(j) + 5) = std50cm;
    end
    clear monthtest;
    % Create quarterly soil temp tests
    ONDtest = (wydatevec_s(:, 2)==10) | (wydatevec_s(:, 2)==11) | ...
        (wydatevec_s(:, 2)==12);
    JFMtest = (wydatevec_s(:, 2)==1) | (wydatevec_s(:, 2)==2) | ...
        (wydatevec_s(:, 2)==3);
    AMJtest = (wydatevec_s(:, 2)==4) | (wydatevec_s(:, 2)==5) | ...
        (wydatevec_s(:, 2)==6);
    JAStest = (wydatevec_s(:, 2)==7) | (wydatevec_s(:, 2)==8) | ...
        (wydatevec_s(:, 2)==9);
    
    quarterscell = {ONDtest, JFMtest, AMJtest, JAStest};
    colindex6 = [75, 81, 87, 93];
    % Then calculate them and  put in soiltempsummary
    for j = 1:4;
        mean5cm = nanmean(nancheck(st5(quarterscell{j})));
        std5cm = nanstd(nancheck(st5(quarterscell{j})));
        mean20cm = nanmean(nancheck(st20(quarterscell{j})));
        std20cm = nanstd(nancheck(st20(quarterscell{j})));
        mean50cm = nanmean(nancheck(st50(quarterscell{j})));
        std50cm = nanstd(nancheck(st50(quarterscell{j})));
        soiltempsummary(i, colindex6(j)) = mean5cm;
        soiltempsummary(i, colindex6(j) + 1) = std5cm;
        soiltempsummary(i, colindex6(j) + 2) = mean20cm;
        soiltempsummary(i, colindex6(j) + 3) = std20cm;
        soiltempsummary(i, colindex6(j) + 4) = mean50cm;
        soiltempsummary(i, colindex6(j) + 5) = std50cm;
    end
    % Calculate some seasonal/yearly means
    soiltempsummary(i, 99) = nanmean(nancheck(st5)); % mast5cm
    soiltempsummary(i, 100) = nanstd(nancheck(st5)); % stdev
    soiltempsummary(i, 101) = nanmean(nancheck(st20)); % mast20cm
    soiltempsummary(i, 102) = nanstd(nancheck(st20)); % stdev
    soiltempsummary(i, 103) = nanmean(nancheck(st50)); % mast50cm
    soiltempsummary(i, 104) = nanstd(nancheck(st50)); % stdev
    % Snowcovered soil temp means
    soiltempsummary(i, 105) = nanmean(nancheck(st5(wysnowtest)));
    soiltempsummary(i, 106) = nanstd(nancheck(st5(wysnowtest)));
    soiltempsummary(i, 107) = nanmean(nancheck(st20(wysnowtest)));
    soiltempsummary(i, 108) = nanstd(nancheck(st20(wysnowtest)));
    soiltempsummary(i, 109) = nanmean(nancheck(st50(wysnowtest)));
    soiltempsummary(i, 110) = nanstd(nancheck(st50(wysnowtest)));
    % Snowfree soil temp means
    soiltempsummary(i, 111) = nanmean(nancheck(st5(~wysnowtest)));
    soiltempsummary(i, 112) = nanstd(nancheck(st5(~wysnowtest)));
    soiltempsummary(i, 113) = nanmean(nancheck(st20(~wysnowtest)));
    soiltempsummary(i, 114) = nanstd(nancheck(st20(~wysnowtest)));
    soiltempsummary(i, 115) = nanmean(nancheck(st50(~wysnowtest)));
    soiltempsummary(i, 116) = nanstd(nancheck(st50(~wysnowtest)));
    
    % Get soil temperatures at transition times (see above)
    if isnan(onsetdatenum_s)
        preonset5cmST = nan;
        preonset5cmSTsd = nan;
        preonset20cmST = nan;
        preonset20cmSTsd = nan;
        preonset50cmST = nan;
        preonset50cmSTsd = nan;
    else
        preonset5cmST = nanmean(nancheck(data{stcol(1)}(preonsettest)));
        preonset5cmSTsd = nanstd(nancheck(data{stcol(1)}(preonsettest)));
        preonset20cmST = nanmean(nancheck(data{stcol(2)}(preonsettest)));
        preonset20cmSTsd = nanstd(nancheck(data{stcol(2)}(preonsettest)));
        preonset50cmST = nanmean(nancheck(data{stcol(3)}(preonsettest)));
        preonset50cmSTsd = nanstd(nancheck(data{stcol(3)}(preonsettest)));
    end
    % Same operations for SoilT before and after melt
    if isnan(meltdatenum_s)
        premelt5cmST = nan;
        premelt5cmSTsd = nan;
        postmelt5cmST = nan;
        postmelt5cmSTsd = nan;
    else 
        premelttest = (datenum_s<=meltdatenum_s & ...
            datenum_s>=(meltdatenum_s-14)); % index the premelt interval
        postmelttest = (datenum_s>=meltdatenum_s & ...
            datenum_s<=(meltdatenum_s+14)); % index the postmelt interval
        premelt5cmST = nanmean(nancheck(data{stcol(1)}(premelttest)));
        premelt5cmSTsd = nanstd(nancheck(data{stcol(1)}(premelttest)));
        postmelt5cmST = nanmean(nancheck(data{stcol(1)}(postmelttest)));
        postmelt5cmSTsd = nanstd(nancheck(data{stcol(1)}(postmelttest)));
        premelt20cmST = nanmean(nancheck(data{stcol(2)}(premelttest)));
        premelt20cmSTsd = nanstd(nancheck(data{stcol(2)}(premelttest)));
        postmelt20cmST = nanmean(nancheck(data{stcol(2)}(postmelttest)));
        postmelt20cmSTsd = nanstd(nancheck(data{stcol(2)}(postmelttest)));
        premelt50cmST = nanmean(nancheck(data{stcol(3)}(premelttest)));
        premelt50cmSTsd = nanstd(nancheck(data{stcol(3)}(premelttest)));
        postmelt50cmST = nanmean(nancheck(data{stcol(3)}(postmelttest)));
        postmelt50cmSTsd = nanstd(nancheck(data{stcol(3)}(postmelttest)));
    end
    % Put in the matrix
    soiltempsummary(i, 117) = preonset5cmST;
    soiltempsummary(i, 118) = preonset5cmSTsd;
    soiltempsummary(i, 119) = preonset20cmST;
    soiltempsummary(i, 120) = preonset20cmSTsd;
    soiltempsummary(i, 121) = preonset50cmST;
    soiltempsummary(i, 122) = preonset50cmSTsd;
    soiltempsummary(i, 123) = premelt5cmST;
    soiltempsummary(i, 124) = premelt5cmSTsd;
    soiltempsummary(i, 125) = postmelt5cmST;
    soiltempsummary(i, 126) = postmelt5cmSTsd;
    soiltempsummary(i, 127) = premelt20cmST;
    soiltempsummary(i, 128) = premelt20cmSTsd;
    soiltempsummary(i, 129) = postmelt20cmST;
    soiltempsummary(i, 130) = postmelt20cmSTsd;
    soiltempsummary(i, 131) = premelt50cmST;
    soiltempsummary(i, 132) = premelt50cmSTsd;
    soiltempsummary(i, 133) = postmelt50cmST;
    soiltempsummary(i, 134) = postmelt50cmSTsd;
end

% Write the soil data files
soilwaterfn = ['wyear_soilwatersummary_' soilinput normstr '.txt'];
soiltempfn = ['wyear_soiltempsummary_' soilinput '.txt'];
csvwrite([processeddatapath soilwaterfn], soilwatersummary);
csvwrite([processeddatapath soiltempfn], soiltempsummary);

junk = 99;

