function checked = nancheck(series)
% nancheck.m
%
% This function checks a timeseries for the presence of NaN's. If more than
% 5% of the input series is NaN's, this returns true (1) so that the
% calling program can determine whether or not to use this data.

lenSeries = length(series);
numNans = sum(isnan(series));
if numNans > 0.05*lenSeries
    checked = nan * zeros(size(series));
else
    checked = series;
end