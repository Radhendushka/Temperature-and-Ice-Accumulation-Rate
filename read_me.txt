The R codes and the data files used in Srivastava and Sengupta (2025), "A semi-parametric model for assessing the effect of temperature on ice accumulation rate from Antarctic ice core data" is given in this repository. 

The two data sets are publicly available from the sites mentioned in the references. 

Vostok Data (Depth, corrected Ice age (GT4) , deltaTS (temperature) ) is available at
https://www.ncei.noaa.gov/pub/data/paleo/icecore/antarctica/vostok/deutnat.txt
Vostok Data (depth, ice age (AICC2012) ) is available at
https://www.ncei.noaa.gov/pub/data/paleo/icecore/antarctica/aicc2012icecore-data.xls

EPICA Dome C data (ztop (depth), Age (EDC3), Temperature) is available at 
https://www.ncei.noaa.gov/pub/data/paleo/icecore/antarctica/epica_domec/edc3deuttemp2007.txt
EPICA Dome C data (depth, ice age (AICC2012) ) is available at
https://www.ncei.noaa.gov/pub/data/paleo/icecore/antarctica/aicc2012icecore-data.xls


Data files

EDC_AICC_model2.csv, Vostok_AICC_model2.csv : This file contains log(AAR) (y), temperature (x), age(z) and depth(d) for locations EPICA Dome C and Lake Vostok. These columns are derived from the downloaded data.



The R script All_functions.R consists of following R functions.  
LS_Estimate: This function computes $\hat{gamma}$ and $\log \hat{g}$ given in section 2. 
Smooth_Estimate: This function computes $\tilde{gamma}$ and $\log {\tilde{g}}$ given in section 2.
BS_Estimate: This function computes the bootstrap standard error, bootstarp quantiles of  $\tilde{gamma}$ and $\log {\tilde{g}}$ given in section 4.  



The R script files Parameter_estimate_EPICA_Dome_C.R and Parameter_estimate_Vostok.R uses the respective data files of the locations and estimates the model parameters along with generate the plots shown in Figure 6, 7 and 8 of the article. The runtime for these scripts are around 75 and 35 minutes respectively on a standard laptop of Intel(R) Core(TM) Ultra 7 155H (1.40 GHz) processor and 32.0 GB RAM.  
