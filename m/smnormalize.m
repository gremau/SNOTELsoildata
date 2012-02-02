function norm=smnormalize(in, normtype, varargin)
% normalization functions for soil moisture timeseries
%
% args:
% in = input 1d array of sensor values
% normtype = integer 1-3 representing normalization type
% varargin = cellarray containing the indexvalue  of in to set as minimum 
%            in normtype 3
%
% Greg Maurer - 2/1/2012

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
    
% normalize by setting min to value at specific time (index given as arg3) 
elseif normtype == 3
    in_min = in(varargin{1});
    norm = (in - in_min);

else
    error('invalid normalization parameter')

end
