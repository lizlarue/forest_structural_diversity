### 1. Run download_cover_data.R script: 
Downloads data for all date site combinations according to the *NEON_sites_dates_for_cover.csv* file from [data folder](../data). Results in *prelim_cover.csv* inside *prelim_cover.rar* compressed file in [output folder](./output). 

### 2. Run make_plot_data_table.R script:
Runs for all date site combinations using *NEON_sites_dates_for_cover.csv* and *All_NEON_TOS_Plot_Centroids_V8.csv* from [data folder](../data). Creates a csv that has plot level descriptions. Uses *helper.R* for coordinate conversion. Results in *plot_data_table.csv* in [output folder](./output). 

### 3. Run download_lidar_data.R script:
Downloads LiDAR data from NEON for all rows in the *plot_data_table.csv* created above, merges all the tiles for each date-site combination, and extract 200\*200m plot files for each easting northing values at those sites. Results in plot level *.laz* files. Not uploaded here, too big files (2.63 GB in total). 

### 4. Run lidar_structural_metrics_extraction.R script:
Runs on each of the plot-level *.laz* files generated above and extracts structural metrics for each one. Results in *lidar_structural_metrics.csv* in [output folder](./output).

### 5. Run plant_cover_metrics_extraction.R script:
Uses the *NEON_sites_dates_for_cover.csv* and *field-sites.csv* from [data folder](../data) and the *plot_data_table.csv* created after Step 2. Extracts plant cover data on a site and plot level and stores it in *cover_by_site.csv* and *cover_by_plot.csv* files respectively in [output folder](./output). <br><br>
Finally, combines the plot-level plant cover and LiDAR data generated in the above step to make *structural_metrics_by_plot.csv* file in [output folder](./output). 

### 6. Run download_aop_data.R script:
Uses *plot_data_table.csv*, and for each plot, downloads hyperspectral data from NEON. Handles cases for single or multiple downloaded files for coordinate mapping to get the reflactance data. Processes the reflectance data and stores the final reflectance values for a plot with buffer in a .h5 format file along with relevant data for further PCA analysis. The final dataset is in AWS S3 bucket in: mukundkalantri.forestbiodiversity/AOP_Data/plot_level_AOP/. Dataset is too big to be uploaded here.

----------------------------

Case to be handled in the files above:
In the download_lidar_data.R and the download_aop_data.R files, handling for the site BLAN needs to be done in a try catch block. This site has data in different UTM zone. The original UTM coordinates needs to be converted to new coordinates based on this UTM zone shift using some R tool/function. Once that is done, all the scripts can be run in the sequence as mentioned above and that should automatically take care of the new files to be generated for this site.
