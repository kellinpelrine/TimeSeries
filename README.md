# TimeSeries

SVAR: structural VAR model examining whether changes in monetary policy have real effects. Loads data from SVAR_Data.csv. Estimates impulse response functions from a monetary policy shock, along with related techniques such as a bootstrap confidence interval. The below figure is one example of output, illustrating the change in GDP (over 50 periods) from a shock to monetary policy.

![image](https://github.com/kellinpelrine/TimeSeries/blob/master/IRF%20with%20Bootstrap%20CI.jpg)



<br/><br/>
Detrend_Corr: processes data from NIPA Table 1.1.3 (from the Bureau of Economic Analysis; .csv included here). Detrends using Hodrick-Prescott filter, examines correlations.

Kendall_Spearman: comparison of sensitivities of Kendall and Spearman stationarity tests. 

Filters: applying some filters to a time series. Produces the below figure:

![image](https://github.com/kellinpelrine/TimeSeries/blob/master/Filters.jpg)

<br/><br/>

ACF_PACF: estimation of ACF and PACF, along with confidence interval. Produces output like the below figure:

![image](https://github.com/kellinpelrine/TimeSeries/blob/master/ACF_PACF.jpg)
