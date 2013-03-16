## binseries.m
A function for binning 2 data series and then generating a means for each bin.

Usage

    [binMean1, binMean2] = binseries(x, y1, y2, topEdge, botEdge, numBins)

* x, y1, y2 = xdata, ydata and ydata2 (usually a std of y1) series
* topEdge = upper x range of bins
* botEdge = lower x range of bins
* numBins = number of bins
* binMean1 = mean of ydata bins
* binMean2 = mean of ydata2 bins

## plot\_examplesites.m
 This makes a number of figures that highlight snow-soil interactions
 at particular sites.
 
 * **File input:** The first 4 figures directly load hourly or daily data. 
 The last 2 use the summary files
 * **User input:** User can select different sites to plot in figs 2, 4,
 5, 6 in the code
 
 Figures:
 
 1. 2 year SWE/temp comparison at Mosby Mtn.
 2. Multi-year SWE/Tsoil/VWC time series at a site (changeable)
 3. Full VWC timeseries of 4 sites - the data for Fig 4
 4. VWC Frequency histograms for 4 contrasting sites
 5. Scatterplot/Regression of snowcov. Tsoil vs pre-onset T at 1 site
 6. Scatterplot/Regression of summer (JAS) VWC vs meltday at 1 site
