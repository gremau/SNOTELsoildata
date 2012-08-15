function m = load7100Avg(desiredAvg)
% load7100Avg.m
%
% Loads long-term average SWE or Precip data from a processed datafile
% and returns the array.
%
% arguments:
% desiredAvg = 'swe' or 'precip'

% ------------------------------------------------------------------------
% Set data path and file name, read in file
datapath = '~/data/current/SNOTELsoil-climate/data_analysis/processed_data/';

if strcmpi(desiredAvg, 'swe')
    datafile = '7100Avg_SWEmm.csv';
    formatstr = ['%u%u%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f' ...
    '%f%f%f%f%f'];
elseif strcmpi(desiredAvg, 'precip')
    datafile = '7100Avg_Precipmm.csv';
    formatstr = '%u%u%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
end

% Generate list site 30yr average data from 7100_avgprecipswe file
fid = fopen([datapath datafile]);
datacell = textscan(fid, formatstr,'Headerlines', 1, 'Delimiter', ',', ...
    'TreatAsEmpty', '*');
fclose(fid);

% Concatenate siteID and all data into m (ID's converted to double)
m = horzcat(double(datacell{2}), datacell{6:end});
% Remove unwanted states and bad data rows (marked with 1 in column 1)
statetest = strcmpi(datacell{4}, 'OR') | strcmpi(datacell{4}, 'WA');
badrowtest = datacell{1}==1;
m((statetest|badrowtest), :) = [];


%junk = 99;
