### 1. Run download_cover_data.R script: 
Downloads data for all date site combinations according to the *NEON_sites_dates_for_cover.csv* file in *forest_structural_diversity/data/*. Results in *prelim_cover.csv* in *prelim_cover.rar* compressed file. \n

### 2. Run make_plot_data_table.R script:
Runs for all date site combinations from *NEON_sites_dates_for_cover.csv* and using *All_NEON_TOS_Plot_Centroids_V8.csv* in *forest_structural_diversity/data/* creates a csv that has plot level descriptions. Results in *plot_data_table.csv*. \n

### 3. Run download lidar_data.R script:
Downloads LiDAR data from NEON for all rows in the *plot_data_table.csv*, merges all the tiles for each date-site combination, and extract 200\*200m plot files for each easting northing values at those sites. Results in plot level *.laz* files. Not uploaded here, too big files (2.63 GB in total). \n

### 4. Run lidar_structural_metrics_extraction.R script:
Runs on each of the plot-level *.laz* files generated above and extracts structural metrics for each one. Results in *lidar_structural_metrics.csv*.
