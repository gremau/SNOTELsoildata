function filled_series = interpseries(series)
% INTERPOLATION to fill gaps in a timeseries. It should fill all NaN's with
% an interpolated value. This can be usefull for running mean and
% variance calculations.
%
% Interpolates over NaNs (data gaps) in the input time series (may be
% complex), but ignores trailing and leading NaN.
%
% WARNING: If series contains only NaN's, interp1 throws an error.
%
% Adapted from FIXGAPS routine on Matlab Central file exchange,
% by R. Pawlowicz 6/Nov/99

if sum(isnan(series))==length(series)
    filled_series = series;
else
    filled_series = series;
    bad = isnan(series);
    good = find(~bad);
    
    bad([1:(min(good)-1) (max(good)+1):end]) = 0;
    
    filled_series(bad)=interp1(good, series(good), find(bad), 'pchip');
end