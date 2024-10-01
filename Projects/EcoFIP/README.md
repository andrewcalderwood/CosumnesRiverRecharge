# README for EcoFIP analyses


# Pre-processing scripts
- AEM_review loads the statewide AEM data and crops to a buffer of the groundwater model extent
- EcoFIP_gw_prep uses the EcoFIP analysis grids to summarize groundwater elevation data to the grid scale
- groundwater_data_review uses the periodic groundwater level measurement database to plot groundwater level hydrographs and contours for the historic period in the study area
- 

# Geologic analysis
- AEM_tprogs_comparison looks at the spatial location of sand/gravel (coarse facies) between the 100 TPROGs realizations and the AEM data on the modflow grid scale
- AEM_coarse_deposits visualizes the the AEM data in the modflow domain as coarse fraction and converts the data to hydraulic conductivity to estimate the potential recharge following the results of Maples et al. 2020
- TPROGS_geomK similarly uses the TPROGS geologic models and depth to water data to estimate recharge potential following Maples et al. 2020
- T2Par Cosumnes loads and cleans the well completion report data to apply it in the Texture2Par code to develop a geologic model that incorporates the AEM data. A test model with the WCR data as percent coarse does not work as it fails to interpolate coarse and fine regions but rather shows very local hot spots.

# Project analysis
- Seepage_analysis summarizes the .sfr.out files from modflow runs for the 10 geologic realizations. It saves the daily and monthly average output to hdf5 files
- seepage_mapping loads the pre-processed hdf5 files to map them onto the EcoFIP grid polygons
- EcoFIP_MF is the start of the script to load the existing modflow model to update the model where the results from the HEC-RAS model will inform the SFR package

# Extra
- model_data_summary was started with the standard flopy model load but I never used it to do anything
