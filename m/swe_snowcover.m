function l = swe_snowcover(site_id, decday_h)
% swe_snowcover.m
%
% Takes swe array and returns a logical matrix indicating snowcover
% for the site
% Data for all years (available) is represented.
%
% ver 1: 110803 GM

%close all;      % clear any figures
%fignum = 0;     % used to increment figure number for plots

% load daily data from site w/ loadsnotel:
[dailyData, ~] = loadsnotel('daily', site_id);

% UNCOMMENT if hourlyData files are not there yet (and comment line above)
% t = [1;1;1;1];
% dailydata = {t '2011-03-29' t t };

% parse out the date/times
decday_d = datenum(dailyData{2}, 'yyyy-mm-dd');

% Snow water equivalent
wteq = dailyData{4}; 

% Convert the daily swe data to an hourly value by matching values with 
% decday_h days and copying to a new array
wteqHourly = zeros(length(decday_h), 1);
for i = 1:length(wteqHourly)
    index = decday_d == floor(decday_h(i));
    wteqHourly(i) = wteq(index);
end

% Logical test for snowcover based on SWE values
swe_snowcover = wteqHourly > 0.5;

l = swe_snowcover;