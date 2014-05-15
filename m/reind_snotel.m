function m_reind = reind_snotel(in, type)
% normalization functions for soil moisture timeseries
%
% args:
% in = input 1d array of sensor values
% normtype = integer 1-3 representing normalization type
% varargin = cellarray containing the indexvalue  of in to set as minimum 
%            in normtype 3
%
% Greg Maurer - 2/1/2012
if strcmpi(type, 'daily')
    tvec = datevec(in{2}, 'yyyy-mm-dd');
    % Create a full daily time vector for the water year
    minyear = min(in{length(in)}); maxyear = max(in{length(in)});
    tstart = datenum(minyear-1,10,1);
    tend = datenum(maxyear,9,30);
    fulltvec = datevec(tstart:tend);
elseif strcmpi(type, 'hourly')
    tvec = datevec(strcat(in{2}, in{3}), 'yyyy-mm-ddHH:MM');
    % Create a full hourly time vector for the water year
    minyear = min(in{length(in)}); maxyear = max(in{length(in)});
    tstart = datenum(minyear-1,10,1);
    tend = datenum(maxyear,9,30,23,0,0);
    fulltvec = datevec(tstart:(1/24):tend);
end

if length(tvec) ~= length(fulltvec)
    
% Get the location in fulltvec where each timestamp of tvec matches
[~, loc] = ismember(tvec, fulltvec, 'rows');
site = unique(in{1});
m_reind{1} = nan * zeros(length(fulltvec),1);
m_reind{1}(:,1) = site;
m_reind{2} = datestr(fulltvec, 'yyyy-mm-dd');
m_reind{3} = datestr(fulltvec, 'HH:MM');

for i=4:length(in)
    m_reind{i} = nan*zeros(length(fulltvec),1);
    m_reind{i}(loc) = in{i};
end

    % Fix the water-year vector at end of matrix
    wyearvec = fulltvec(:, 1);
    wytest = (fulltvec(:,2)==10 | fulltvec(:,2)==11 | fulltvec(:,2)==12);
    wyearvec(wytest) = wyearvec(wytest) + 1;
    m_reind{length(in)} = wyearvec;

else
    m_reind = in;

end

