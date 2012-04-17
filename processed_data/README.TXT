SNOTELsoil-climate/data_analysis/processed_data

This directory contains files that were created by data analysis scripts or
spreadsheets. From here they are read by other scripts and the data is used in
plotting or for other uses. 

*** BE VERY CAREFUL changing these files!!!***

Any changes to files themselves (such as column additions or deletions) is
likely to break the data analysis process that use them. Only update these
files by changing the source script/spreadsheet and then re-running/exporting
document changes, and make necessary changes to the scripts that rely on these 
files after they are changed.

ORIGIN OF THESE FILES

7100_*.csv : 
Exported from rawdata/longterm_averages/SnowSurvey_7100Avg_MASTER.gnumeric
Elevation data in these files are unreliable

SNOTELinventory.csv :
Exported from rawdata/station_inventory/NRCS_inventory_MASTER.gnumeric

SNOTELrangelist.csv :
Not sure where this came from, but it is used in at least one script

wyear_*.csv :
Generated with m/summarize_wateryears.m script
