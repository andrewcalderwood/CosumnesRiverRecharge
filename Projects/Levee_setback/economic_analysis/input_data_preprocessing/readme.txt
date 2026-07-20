# README for input preprocessing

00_dem_soil_preparation and 00_hydrologic_data jupyter notebook are both files that are only needed once except if the model is being extended into the future such as 2025-2026 then hydrologic data may need to be re-run.

01_rep_SWB_runs is the script that creates the optimized SWB for a pre-defined range of DTW values, 0 to 300 ft, and for cases with and without POD by each crop. These then create look up tables that will speed up the connected model runs as they will only need to interpolate results. This would break if in the future projects, there are DTW> 300 ft, but this would also be problematic as many wells are in the 200-300 ft range.

01_rep_SWB_runs-Copy1 is just a temp verison to re-run the first year for the baseline that got overwritten by accident.

These files should not need to be re-run unless an error in the SWB optimization is found or if there is a need to update DTW range or other input parameters such as the climate range.

