% ## binseries.m
% A function for binning 2 data series and then generating a means for each bin.
%
% Usage
%
%    [binMean1, binMean2] = binseries(x, y1, y2, topEdge, botEdge, numBins)
%
% * x, y1, y2 = xdata, ydata and ydata2 (usually a std of y1) series
% * topEdge = upper x range of bins
% * botEdge = lower x range of bins
% * numBins = number of bins
% * binMean1 = mean of ydata bins
% * binMean2 = mean of ydata2 bins

function [binMean1, binMean2] = binseries(x, y1, y2, topEdge, botEdge, numBins)
    addpath('/home/greg/data/code_resources/m_common/nanstats/');
    binEdges = linspace(botEdge, topEdge, numBins+1);
    [h,whichBin] = histc(x, binEdges);
    for i = 1:numBins
        flagBinMembers = (whichBin == i);
        binMembers1     = y1(flagBinMembers);
        binMembers2     = y2(flagBinMembers);
        binMean1(i)     = nanmean(binMembers1);
        %binStd1(i)      = nanstd(binMembers1);
        binMean2(i)     = nanmean(binMembers2);
        %binStd2(i)      = nanstd(binMembers2);
    end
end
