***About***
The scripts in this folder (using MATLAB 2016a, Python 2.7.10, Jupiter notebooks 4.1.0, Landsat-util 0.12.1, and GDAL framework package 1.11) are used to acquire Landsat imagery, manage files, (calculate ice velocity using imager cross-correlation,) reproject data, filter data, and produce velocity mosaics and velocity time series. The code here is tailored to studying the Larsen C Ice Shelf on the western Antarctic Peninsula.

***Citation***
Pope, Allen (2016). Processing Landsat 8 Velocities for Larsen C. Zenodo. _doi________
Available at: https://github.com/allenpope/Landsat8_Velocity_LarsenC

***Workflow Information***
Scripts are designed to be used in the following order:
	1) Download_Prep_L8.ipynb
	2) PyCorr (see Fahnestock et al.)
	3) Subset_pycorr.ipynb
	4) Larsen_to_nc.m
	5) Larsen_filter_nc.m
	(Iteration based on one scenes, replacing "93Larsen_filter_nc.m" with "Larsen_update_nc.m")
	6) Larsen_mosaic.m / Larsen_timeseries.m

***Literature Citations***
Borstad, C., McGrath, D., & Pope, A. (2016). Rapid growth of a large rift on the Larsen C Ice Shelf: Observations and Implications. Geophysical Research Letters, in review.

Fahnestock, M, Scambos, T., Moon, T., Gardner, A., Haran, T., & Klinger, M. (2016). Rapid large-area mapping of ice flow using Landsat 8. Remote Sensing of Environment, 185, 84-94. http://dx.doi.org/10.1016/j.rse.2015.11.023
