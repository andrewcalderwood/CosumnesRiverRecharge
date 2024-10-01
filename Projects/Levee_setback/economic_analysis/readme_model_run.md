# Model Run README

# Pre-processing
Pre-processing is done to inform set up of the soil water budget, irrigation optimization, and crop choice models. Additionally there is pre-processing needed to translate results from these models to the MODFLOW model and back.

The scripts that support this include: 
- 00_dem_soil_preparation
- 00_hydrologic_data

# Model testing
Several scripts were developed in the model testing phase to create the soil water budget and the irrigation optimization version of it

The scripts for this include:
- 01_SWB_testing to run the optimization for a crop and year of interest over multiple fields to test run times/optimizer settings
- 02_iterative_irrigation_optimizer which was set up to run the irrigation optimizer over multiple crops for a specific year of interest to save the output for each year of output for percolation and irrigation (gw pumping). The batch file run_irrigation_optimizer.bat was used to help with testing of this.

# Model Running
The final connected model was set up with a variety of pre-processed inputs and functions to go from a start year to an end year by using the existing modflow model and re-writing the RCH/WEL inputs using the output from the crop choice to inform the optimized soil water budget plus two more soil water budgets for native and ag fields during winter without optimization.

03_model_connect references the following scripts generally in this order:
- from functions.f_gw_dtw_extract import sample_dtw, avg_heads to sample groundwater depth to water from the MODFLOW output of the previous step
- parcelchoicemodelupdate.f_predict_landuse for the crop choice model
- functions.Basic_soil_budget_monthly for soil water budget functions
- from f_rep_swb_profit_opt import load_run_swb to run the optimizer over the soil water budget for representative depth to water profiles
- from functions.output_processing import get_local_data, out_arr_to_long_df, get_wb_by_parcel to load the output from the optimized soil water budget to input to RCH/WEL
- from reference_swb_ag_winter import run_swb_ag_winter to run the soil water budget model in the winter to estimate deep percolation for ag fields with the crops for the current year

# Post-processing
03b_summarize_output takes the saved irrigation rates and the simulated DTW after the fact to calculate the actual profit and yield on a parcel basis. This code could be integrated into 03_model_connect if interested in getting year by year updates in profit to inform decision making. This would also help provide better estimates of percolation.




