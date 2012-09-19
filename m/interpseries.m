function filled_array = interpseries(array)
% INTERPOLATION to fill gaps in a timeseries. It should fill all NaN's with
% an interpolated value. This can be usefull for running mean and
% variance calculations.
%
% Interpolates over NaNs (data gaps) in the input time series (may be
% complex), but ignores trailing and leading NaN. 
% WARNING: If array contains only NaN's, interp1 throws an error.
%
% Adapted from FIXGAPS routine on Matlab Central file exchange,
% by R. Pawlowicz 6/Nov/99

filled_array = array;

bad = isnan(array);
good = find(~bad);

bad([1:(min(good)-1) (max(good)+1):end]) = 0;

filled_array(bad)=interp1(good, array(good), find(bad), 'pchip');