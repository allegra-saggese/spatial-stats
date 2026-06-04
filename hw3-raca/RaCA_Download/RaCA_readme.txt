####################################################################################
List of files and subfolders of RaCA_Download
Download and unzip RaCA_Download.zip 
 
Available on Cloudvault:
https://www.cloudvault.usda.gov/index.php/s/nZDwWAdlmKUW1g4

In the RaCA_Download folder - subfolders and files used for 
NRCS RaCA Summary Reports

To USE: in R script change the  source to the location of the 
unzipped RaCA_Download.zip or individual files

###################################################################################
Background - folder with background documents about how RaCA data was designed, 
             collected and measured.

	Alter_site_location.pdf
	Location Insructions.pdf
	RaCA Data and Field Collection Protocols.pdf
	Rapid Carbon Assessment Field Laboratory Instructions.pdf
	Rapid Carbon Assessment Workbook 2 15.xlsx
	Region_soil components GROUPS.xlsx


############################
Data - Data files (csv) and R scripts include sample data 
	(Bulk Density and Carbon Concentraion) 
	files used in manipulation and interim data files to ease reproduction

	RaCA_data_columns.csv    - explanation of RaCA_sample and RaCA pedon columns
	RaCA_soilDB_elements.csv - explanation of elements available through 
				   fetchRaCA() in soilDB
	RaCA_samples.csv         - raw RaCA data, all samples including lab analysis
                                   and bulk density modeled values
	LUGR_pixelcount.csv      - created from a grid (30m x 30m) of NLCD and 
                                   gSSURGO products, clipped by RaCA region

	Data\BD_model.zip
		RaCA_new_bd_final.R       - will recreate modeled bulk density and 
                                            graphs comparing measured and modeled 
                                            bulk density.
		RaCA_samples_July2016.csv - raw RaCA data, all samples
		NCSS_SOC_abovebelow.csv   - data derived from the NCSS 
                                            Characterization Dataset
		new_horizon.csv           - General Horizon Table used to  
                                            standardizeand relabel 
                                            horizon designations
		RF_O_all.RData (random forest model - predicts bulk density 
                                for organic horizons)
		RF_sab.RData  (random forest model - predicts bulk density 
                                for mineral horizons)
	
	      output: RaCA_samples.csv, BD_resid_eval.csv


	Data\Summary
		RaCA_SOCstock_calc.R - script uses sample and BD model output 
                                       to calculate pedon stocks
	     
	      output: RaCA_SOC_pedons.csv
	
		RaCA_SOC_LUGR&SUmmary.R  - script uses pedon output and pixel counts  
                                           by LUGR classes to calculate weighted 
                                           averages by various classes.
		
	      output:  	RaCA_LUGR_stocks.csv - average value of each LUGR
			bt_RaCAall_stocks.csv - log transformed, averaged then back 
                                                transformed LUGR values,
                                                filled with LU and Group means
					        product used to map SOC stocks 
                                                across CONUS 
		 	All_wt_avg.csv - weighted averages by overall, MO and LULC  
                                         class(all stacked in one file)
			wt_avg/MO_wm_all_wSD.csv - weighted averages for each region
								

	Data\Permutations
		lugr_prob.R - calculated selection probabilities of each pedon 
                              (based on pixel counts of each mapped LUGR class)
		
		SOC_stocks_perm_errors_w_triangle.R - permutations of measurement 
                                                      errors for each sample; 
                                                      summarized by pedon, 
                                                      then LUGR and other classes
		Horizon_propSD.csv
		BD_resid_eval.csv
		LUGR_prob.csv
		

	      output: Data\Permutations\Perm_output
			mc_5_stocks.csv
			mc_30_stocks.csv
			mc_100_stocks.csv
				


			
##################################################################################
Remainder of work was done in ArCGIS

To create smoothed CONUS map:
#work done for first RaCA project#
- Began with NLCD 2011 snapped to gSSUGO 2012 (30m)
- MO boundary 2008 was used to clip each gSSURGO to 17 CONUS MOs 
   (also referred to as RaCA regions)
- Each grid value was decomposed to NLCD class and MUKEY 
- RaCA LUGR classes were assigned to each (based on NRI LULC and RaCA Soil Groups)
- Geometric means (t_LUGRall_wCL.csv) of LUGR classes, including fills of 
    overall LU and MO, b

#new RaCA work#
- The LUGR means were attached to grid cells and mapped
- the SOC pedon stocks were plotted by x y location
- Ordinary Kriging was used to interpolate between points
- The LUGR grids and interpolated mapps (50% transparent) were overlain for 
   presentation
- The MO weighted means (from ALL_wt_avg.csv) were attached to MOs/regions

