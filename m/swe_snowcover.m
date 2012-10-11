function returned = swe_snowcover(site_id, decday_h)
% swe_snowcover.m
%
% Takes a site an hourly decimal day array array and returns a boolean
% matrix indicating HOURLY snowcover for the site
% Data for all years (available) is represented.
%
% ver 1: 110803 GM

% load daily data from site w/ loadsnotel:
dailyData = loadsnotel(site_id, 'daily', 'exclude');

% Round hourly datenum down to days
hourlydays = floor(decday_h);

% Parse out the daily datenum and swe
decday_d = datenum(dailyData{2}, 'yyyy-mm-dd');
wteq = dailyData{4}; 

% Convert the daily swe data to an hourly value by matching values with 
% decday_h days and copying to a new array
wteqHourly = nan * zeros(length(decday_h), 1);
for i = 1:length(decday_d)
    hourlyindex = hourlydays==decday_d(i);
    wteqHourly(hourlyindex) = wteq(i);
end

% Logical test for snowcover based on SWE values
swe_snowcover = wteqHourly > 0.1;

returned = swe_snowcover;
