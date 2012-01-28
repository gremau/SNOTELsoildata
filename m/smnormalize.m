function norm=smnormalize(in, normtype)
% normalization functions for soil moisture timeseries

% 0-1 normalization (0=min and 1=max for series)
if normtype == 1
    in_max = max(in);
    in_min = min(in);
    norm = (in - in_min)./(in_max-in_min);

% subtract min value, all series have a zero point (problem if there are
% negative values)
elseif normtype == 2
    in_min = min(in);
    norm = (in - in_min);
    
% normalize min to value at specific date (day 263.1) for every timeseries 
elseif normtype == 3
    in_min = in(12495);
    norm = (in-in_min);

else
    error('invalid normalization parameter')

end