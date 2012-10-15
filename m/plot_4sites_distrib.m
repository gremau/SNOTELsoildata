% plot_4sites_distrib.m
%
% Plots quarterly soil moisture distributions of 4 SNOTEL sites for all
% (available) years.
%
% Greg Maurer 120202

close all; clear all;   % clear figures and variables in the workspace
fignum = 0;     % used to increment figure number for plots
addpath('/home/greg/data/programming/m_common/');

% Set list of sites, sensor output(vwc or temp), and sensor depth
siteIDs = [828, 333, 452, 336]; % elev/swe: hi/hi, low/hi, hi/low, low/low
% 828 = TrialLk, 333 = BenLomTrail, 452=DonkeyRes, 573=Big BendNV
sensoroutput = 'vwc';
sensordepth = 2; %(1=5cm, 2=20cm, 3=50cm);
startwy = 2006;

% Select TEMP or VWC data and set distribution bins and plot axes
if strcmpi(sensoroutput, 'vwc');
    sensorcolumn = sensordepth + 3; % get proper column using sensordepth
    xmin = 0;
    % If running RAW SENSOR DATA (no normalization)
    % xedges = 0:1:100; % raw sm data bins (0-100)
    % xmax = 75
    % If running NORMALIZED data with smnormalize
    xedges = 0:0.01:1; % normalized vwc bins (0-1)
    xmax = 1; % these axes are good for normalized data
    ymax = 0.125;
    disp('*** Running in normalized soil moisture data mode ***');
elseif strcmpi(sensoroutput, 'temp');
    sensorcolumn = sensordepth + 6; % get proper column using sensordepth
    xedges = -15:0.5:35; % soil temp bins -15-35 deg C
    xmax = 35; % corresponding x axis limits
    xmin = -15;
    ymax = 0.2;
end

% Set up PLOT 1 - add each site's timeseries on iteration through following
% loop
fignum = fignum+1;
h = figure(fignum);

% Allocate for histogram and mean matrices - fill in following loop
histograms = zeros(length(xedges), 16);
means = zeros(16, 1);

for i = 1:length(siteIDs);
    % Load hourly data from site  w/ loadsnotel:
    siteHourly = loadsnotel(siteIDs(i), 'hourly', 'exclude');
    % Get rid of wateryears prior to startwy
    wyexclude = siteHourly{10}>startwy-1;
    for j = 1:10
        siteHourly{j} = siteHourly{j}(wyexclude);
    end
    % Parse out the desired sensor depth, normalize if plotting vwc
    if strcmpi(sensoroutput, 'vwc')
        sensordata = filterseries(siteHourly{sensorcolumn}, 'sigma', 25, 3);
        % SPECIAL Case for Taylor Cyn - There is some bad data that makes it
        % past filter and messes up normalization - remove it
        if siteIDs(i) == 336 %811n for TaylorCyn, 336 for BigBend
            test = sensordata<10; %18 for TaylorCyn, 10 for BigBend
            sensordata(test)=nan;
        end
        sensordata = smnormalize(sensordata, 1);
    else
        sensordata = siteHourly{sensorcolumn};
    end
    
    % Create date arrays
    datevec_h = datevec(strcat(siteHourly{2}, siteHourly{3}), 'yyyy-mm-ddHH:MM');
    datenum_h = datenum(datevec_h);
    
    % PLOT 1 - add the entire timeseries for site i
    ticklocations = linspace(min(datenum_h), max(datenum_h), 20);
    set(h, 'Name', ['Site ' num2str(siteIDs(i)) ' - Full sensor timeseries']);
    subplot (4, 1, i)
    plot(datenum_h, sensordata, 'k');
    set(gca, 'XTick', ticklocations);
    set(gca, 'XTickLabel', ticklocations);
    datetick('x', 12, 'keepticks');
    ylabel(sensoroutput);
    title(['Site ' num2str(siteIDs(i))]);
    
    % Create logical tests and pull desired quarters (3 months intervals)
    testOND = (datevec_h(:,2)==10 | datevec_h(:,2)==11 | datevec_h(:,2)==12);
    testJFM = (datevec_h(:,2)==1 | datevec_h(:,2)==2 | datevec_h(:,2)==3);
    testAMJ = (datevec_h(:,2)==4 | datevec_h(:,2)==5 | datevec_h(:,2)==6);
    testJAS = (datevec_h(:,2)==7 | datevec_h(:,2)==8 | datevec_h(:,2)==9);
    sensordata_OND = sensordata(testOND,:);
    sensordata_JFM = sensordata(testJFM,:);
    sensordata_AMJ = sensordata(testAMJ,:);
    sensordata_JAS = sensordata(testJAS,:);
    
    % Generate histograms for each quarter's sensor data
    histOND = histc(sensordata_OND, xedges);
    histJFM = histc(sensordata_JFM, xedges);
    histAMJ = histc(sensordata_AMJ, xedges);
    histJAS = histc(sensordata_JAS, xedges);
    
    % Normalize and put histograms in the histograms matrix
    % (in plotting order)
    histograms(:, i) = histOND./sum(histOND);
    histograms(:, i+4) = histJFM./sum(histJFM);
    histograms(:, i+8) = histAMJ./sum(histAMJ);
    histograms(:, i+12) = histJAS./sum(histJAS);
    
    % Calculate mean and standard deviations of data from each quarter and
    % place in appropriate vector (in plotting order)
    means(i) = mean(sensordata_OND(~isnan(sensordata_OND)));
    means(i+4) = mean(sensordata_JFM(~isnan(sensordata_JFM)));
    means(i+8) = mean(sensordata_AMJ(~isnan(sensordata_AMJ)));
    means(i+12) = mean(sensordata_JAS(~isnan(sensordata_JAS)));
    % stddevOND = std(sensordata_OND(~isnan(sensordata_OND)));
    % stddevJFM = std(sensordata_JFM(~isnan(sensordata_JFM)));
    % stddevAMJ = std(sensordata_AMJ(~isnan(sensordata_AMJ)));
    % stddevJAS = std(sensordata_JAS(~isnan(sensordata_JAS)));

end
    clear testOND testJFM testAMJ testJAS;


%-------------------------------------------------------------
% PLOT 2. Plot quarterly distributions for all sites
titles = {['Site ' num2str(siteIDs(1)) ' Oct-Dec'] ...
    ['Site ' num2str(siteIDs(2)) ' Oct-Dec']...
    ['Site ' num2str(siteIDs(3)) ' Oct-Dec']...
    ['Site ' num2str(siteIDs(4)) ' Oct-Dec'] ...
    'Jan-Mar' 'Jan-Mar' 'Jan-Mar' 'Jan-Mar' ...
    'Apr-Jun' 'Apr-Jun' 'Apr-Jun' 'Apr-Jun' ...
    'Jul-Sep' 'Jul-Sep' 'Jul-Sep' 'Jul-Sep'};
fignum = fignum+1;
h = figure(fignum);
set(h, 'Name', ['4 SNOTEL Sites - ' sensoroutput ...
    ' quarterly histograms - all years combined']);
% Loop through 16 subplots and plot histograms and means
for i = 1:16;
    subplot (4, 4, i)
    bar (xedges, histograms(:, i), 'g');
    hold on
    plot([means(i) means(i)], [0 1], ':k');
    axis([xmin xmax 0 ymax]);
    ylabel('Frequency');
    title(titles{i});
end