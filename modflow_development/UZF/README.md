This folder contains scripts to pre-process information related to the soil water budget such as evapotranspiration, soil properties and topographic characteristics.  

The scripts are generally run in this order:
1. `Simple_crop_input` loads landuse data, crop coefficient database and potential ET to calculate the crop ET. It originally produced only ETc for each grid cell, and was updated to provide simply the ETc for the list of crops to then input to the soil budget. This script is also used to save a cleaned version of the land use data.  
    - the old version of this file is `UZF_crop_input_and_analysis` which was based on USDA NASS data. It produced gridded ETc for all model cells with different landuse for each year.  
2. `SSURGO processing` translates the SSURGO shapefile and text file to the model grid format
    - `field soil preparation` translates the SSURGO data to the land use polygon format
3. `Basic soil budget` calculates the soil water budget based on IDC equations with simplification:
    - `gridded` has the calculation for each model grid cell based on the gridded ETc and soil data
    - `fields` has the calculation for each agricultural and native land use field
4. `GDE EVT` uses shapefiles of Groundwater Dependent Ecosystem locations created by The Freshwater Trust for the South American Subbasin GSP to identify locations for groundwater evatranspiration in the flow model. 