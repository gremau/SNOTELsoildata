function l = soiltemp_snowcover(Ts)
% soiltemp_snowcover.m
%
% Takes soil temp array and returns a logical matrix indicating snowcover
% for the site
% Data for all years (available) is represented.
%
% ver 1: 110328 GM

%close all;      % clear any figures
%fignum = 0;     % used to increment figure number for plots

sd_threshold = 0.1;
hoursinday = 24;
number = length(Ts);

%
% Calculate a DISCRETE 24hour mean/stdev (requires decday array)
%
% Pad soil temp and date array with nans to an even multiple of 24
% Ts(end + (1:(ceil(number/hoursinday)*hoursinday-number))) = nan; 
% decday(end + (1:(ceil(number/hoursinday)*hoursinday-number))) = nan;
%
% Reshape both into columns of 24 for averaging
% Ts_24 = reshape(Ts, hoursinday,[]);
% decday_24 = reshape(decday, hoursinday,[]);
%
% Logical array of standard deviations
% discstd = std(Ts_24);


%
% Calculate a 24hour RUNNING mean/stddev
%
int = 24;
runmean = filter(ones(int,1)/int,1,Ts);
runvar = (filter(ones(int,1),1,Ts.^2) - int * runmean.^2)/(int-1);
runstd = sqrt(runvar);

% Create logical array of hourly/daily values less than sd_threshold
l = runstd < sd_threshold;

