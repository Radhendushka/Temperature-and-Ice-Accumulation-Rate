The R codes and the data files used in Srivastava and Sengupta (2025), "A semi-parametric model for assessing the effect of temperature on ice accumulation rate from Antarctic ice core data" is given in this repository. 

The two data sets are publicly available from the sites mentioned in the references. 
The AICC2012 age scales for Lake Vostok and EPICA Dome C (aicc2012icecore-data.xls) were downloaded from https://www.ncei.noaa.gov/access/paleo-search/study/15076. 

Data files


EDC_AICC_model2.csv, Vostok_AICC_model2.csv : This file contains log(AAR) (y), temperature (x), age(z) and depth(d) for locations EPICA Dome C and Lake Vostok. These columns are derived from the downloaded data.



The R script file All_functions.R consists of several R functions. 
The function LS_Estimate computes the least square estimates of the model parameters. 
The function Smooth_Estimate computes the proposed smooth estimates of log(g), and corresponding estimates of γ, β and 𝜎.
The function BS_Estimate generate the model based bootstrap samples as described in the article and compute the bootstrap se, 2.5% and 97.5% bootstrap quantiles of the model parameters.  

The R script files rho_n_EPICA_Dome_C.R and rho_n_Vostok.R uses the data files of the respective locations and estimates ϱ along with the generates the plots shown in Figure 6 of the article.


The R script files Parameter_estimate_EPICA_Dome_C.R and Parameter_estimate_Vostok.R uses the respective data files of the locations and estimates the model parameters along with generate the plots shown in Figure 7 and 8 of the article.
