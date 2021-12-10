# Overlap correction algorithm #

Repository for the overlap correction for CHM15k. This code was used for the paper submitted in AMT.


## Content ##
* Calculate overlap correction funtion ( script_overlap_routine_v2.m calling  calculate_overlap_automatic_structured_2.m)
* plot all results (read_overlap_cor_v5). Including:

1. All daily overlap functions
1. temperature model
1. Impact on gradient and temperature model simple PBL detection

* Gradient Analysis (read_gradient)


## Detailed dexcription of standard way of using ##
	1. script_overlap_routine_v3_EPROF.m
		a. fill variables : 
			* stn: wigos-id
			* start_time/end_time: time period for which to calculate overlap correction
			* folder_data : location of E-PROFILE netCDF level 1 files   (NOTE: you'll prefer to transfer ncdf files directly on your machine to improve execution speed (locally ~12hours/1year of data))
			* folder_out :  output directory   (NOTE: prefer a local directory on your machine for speed)
		b. run script                              
	2. read_overlap_cor_v6_eprofile.m
		a. fill variables:
			* station: wigos-id
			* folder_ncdata:  location of E-PROFILE netCDF level 1 files
			* folder_correction:  File path of the output created by script_overlap_routine_v3_EPROF.m in step 2)
			* folder_output: output directory  (NOTE: prefer a local directory for speed)
		b. If you are configuring a new station add a case in the station inputs switch. Then specify:
			* timerange of the correction files (in info.xxx) + timerange of daily correction visualization for actually testing the correction on measured profiles (in info_test.xxx).
			* optical module id (info.tub) and instrument serial number (info.chm). You'll find this information in daily netcdf files attributes.
		c. run script             
		d. filter bad quality overlap estimates
			* identify in plots bad quality dates (best done in Matlab figure by opening >view>plot browser)
			* add these bad days to list_dates_bad_quality in the station switch
		e. run script again


## Examples ##
### Raw ###
![Raw](https://bitbucket.org/repo/zb4zbB/images/1504515377-PR2_grad_20140715_all_raw.png =50x50)
### Corrected with model ###
![Model](https://bitbucket.org/repo/zb4zbB/images/3933804104-PR2_grad_20140715_all_model%20correction.png =50x50)