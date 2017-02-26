To update the data

1. Update the 3-10 yrs nominal and TIPS yields data in worksheet "Merged 1999 - Stata" of "Data_1999.xls" using the following links:
Nominal Yields
http://www.federalreserve.gov/econresdata/researchdata/feds200628.xls
TIPS Yields
http://www.federalreserve.gov/econresdata/researchdata/feds200805.xls

2. Update the 3-month T-Bills yield data in worksheet "3mon Nominal weekly 1999-" of "Data_1999.xls" using the following link:
3-month Nominal Yields
http://www.federalreserve.gov/datadownload/Output.aspx?rel=H15&series=bf876e40d4bdea046af749aae2168603&lastObs=&from=&to=&filetype=csv&label=include&layout=seriescolumn

NOTE: There are weeks when the final weekly observation for the nominal yields is available on a Friday and the final weekly observation for the TIPS yields is available on a Thursday, or vice versa. In those cases, pick the final weekly observations for both and assign them to the Friday of the week.

NOTE: TIPS data for the 3-yr and 4-yr maturities from the beginning of 1999 to the end of 2003 were added from a data file provided by Tobias Adrian, since these data were not available on the board website.



To run the model

3. Install the MFE toolbox available at http://www.kevinsheppard.com/MFE_Toolbox

4. Set the start and end dates for the weekly data in line 20 and 21 of daily2weekly.m

5. Run adrian_wu_rep.m

NOTE: The TIPS market was quite illiquid in the initial years when they were first introduced in the US in 1999. When estimating the model with data going back to 1999 rather than 2013, some parameters in the model that have been restricted to 0 need to be allowed to be different from 0.

NOTE: It took 9596 seconds to estimate the parameter of the model (2003-2010) on a 3.6 GHz 16GB RAM machine.