function filterseries_test()
% filterseries_test.m
%
% This scripts tests the filterseries.m function. Snotel sites are loaded
% and filtered and the results are plotted
%
% ver : 110808 GM

close all;      % clear any figures
fignum = 0;     % used to increment figure number for plots
addpath('/home/greg/data/rawdata/SNOTEL_data/matlab/');

% Set data path and file name, read in file
datapath = '/home/greg/data/rawdata/SNOTEL_data/';

% List of sites to analyze
havedata = [828, 972, 332, 333, 628, 766, 596, 896, 971, 1054, 474, 0];

% Or generate a list from the _sitelist.txt file
%havedata = unique(dlmread([datapath 'soilsensors_hourly/_sitelist.txt']));
% 
% Generate list of sites and their elevations from inventory file
% Create format string (station,elev,cdbs_id only here)
% formatstr = '%*s%f%*s%*s%*s%*s%f%*s%*s%*s%*s%*s%*s%*s%*s%*u';
% fid = fopen([datapath 'station_inventory/UT_soilstations.csv']);
% listcell = textscan(fid, formatstr,'Headerlines', 1, 'Delimiter', ',');
% fclose(fid);
% 
% sitesarray = [listcell{1}, listcell{2}];
% clear listcell;

% Filter out sites for which there is no data in the data folder
% havedatatest = ismember(sitesarray(:, 1), havedata);
% sitesarray = sitesarray(havedatatest, :);
% 
% % Remove bad sites - only 583 so far
% badsitetest = sitesarray(:,1) == 583;
% sitesarray(badsitetest, :) = [];

%elev = [9992, 6700, 8000, 5829, 8967, 9640, 8200, 7250, 7740, 6779, 8000, 8300];
% list = sitesarray(:,1);
% elev = sitesarray(:,2);
%list = [828];
%elev = [9992];
% t = zeros(13, length(list));
% template = {t t t t t t};
% [snowmeans freemeans snowsums freesums snowdev freedev] = template{:};

m = cell(1, length(havedata));

for i = 1:length(havedata);
    [m, ~] = loadsnotel('hourly', havedata(i));
    originalSeries = m{7};
    
    % INTERPOLATION to fill gaps in Ts (helps with running mean and
    % variance calculations below.
    %
    % Interpolates over NaNs (data gaps) in the input time series (may be
    % complex), but ignores trailing and leading NaN.
    %
    % from FIXGAPS routine on Matlab Central file exchange,
    % by R. Pawlowicz 6/Nov/99

    filledSeries = originalSeries;

    bad = isnan(originalSeries);
    good = find(~bad);

    bad([1:(min(good)-1) (max(good)+1):end]) = 0;

    filledSeries(bad)=interp1(good, originalSeries(good), find(bad), 'pchip');
    
    originalSeries = filledSeries;
    %
    % Filter
    %Ts = filterseries(Ts, 'shift', 2.5);
    % FILTER Tsoil data (returns filtered and re-interpolated array) 
    meandiffSeries = filterseries(originalSeries, 'mean', 7);
    %shiftMeandiffSeries = filterseries(Ts_meandiff, 'shift', 5);
    
    % PLOT original data over interpolated data to view differences
    h = figure;
    set(h, 'Name', ['Site ' num2str(havedata(i)) ' - filterseries.m data interpolation']);
    plot(originalSeries, '.r');
    hold on
    %plot(filteredFilledSeries, '.b');
    plot(meandiffSeries, '.k');
    title('Filtered points in red, interpolated data in blue')
end
