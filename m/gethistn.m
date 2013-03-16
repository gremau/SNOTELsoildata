function histn = gethistn(sitelist)

sensordepth = 2; %(1=5cm, 2=20cm, 3=50cm);
startwy = 2006;

% Set distribution bins and plot axes
sensorcolumn = sensordepth + 3; % get proper column using sensordepth
xmin = 0;
% If running RAW SENSOR DATA (no normalization)
% xedges = 0:1:100; % raw sm data bins (0-100)
% xmax = 75
% If running NORMALIZED data with smnormalize
xedges = 0:0.02:1; % normalized vwc bins (0-1)
xmax = 1; % these axes are good for normalized data
ymax = 0.25;

a = []; %92*24*length(sites_hihi);
for i = 1:length(sitelist);
    % Load hourly data from site  w/ loadsnotel:
    siteHourly = loadsnotel(sitelist(i), 'hourly', 'exclude');
    % Parse, filter, and normalize the desired sensor data
    sensordata = filterseries(siteHourly{sensorcolumn}, 'sigma', 25, 3);
    sensordata = smnormalize(sensordata, 1);
    
    % Create date arrays
    datevec_h = datevec(strcat(siteHourly{2},siteHourly{3}),...
        'yyyy-mm-ddHH:MM');
    datenum_h = datenum(datevec_h);
    
    % Get rid of wateryears prior to startwy
    wyexclude = siteHourly{10}>startwy-1;
    sensordata = sensordata(wyexclude);
    datevec_h = datevec_h(wyexclude,:);
    
    % SPECIAL Case for Taylor Cyn - There is some bad data that makes it
    % past filter and messes up normalization - remove it
    %         if siteIDs(i) == 336 %811n for TaylorCyn, 336 for BigBend
    %             test = sensordata<10; %18 for TaylorCyn, 10 for BigBend
    %             sensordata(test)=nan;
    %         end
    
    % Pull desired quarters (3 months intervals) with a logical
    % test
    testJAS = (datevec_h(:,2)>6 & datevec_h(:,2)<10);
    a = [a; sensordata(testJAS)];
end
histogram = histc(a, xedges);
histn = histogram./sum(histogram);
end