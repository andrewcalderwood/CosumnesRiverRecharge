# historical streamflow model development notes

# model results interpretation

right now the model is showing large fluctuations in groundwater level, but I just realized this is because it is still controlled by the specific storage since the model is confined. i could activate the storativity option to estimate that instead to reduce fluctuation. Or by running the model confined I need to increase specific storage

Start by only calibrating Ss because Sy won't have an impact

The calibration could be getting over parameterized at this point. But I may need to add a second scaling factor to the SFR VKA to account for different vertical anisotropy at the stream scale and improve that representation, but this also depends greatly on which TPROGs realization I use. Unless I return to using the sediment data from Constantine and calibrate sand/gravel vs silt/mud. This is going to be the key feature of interest of whether I can accurately represent stream-aquifer interaction. Once I get regional groundwater levels to the appropriate magnitude and fluctuation with recharge and pumping then I need to look at scaling stream bed conductivity using the streamflow data available.