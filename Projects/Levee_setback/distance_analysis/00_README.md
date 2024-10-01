# Order of code processing
1. 08_connec3d creates the input files for Connec3D which must be run then post-processed to identify the connected paths from above land surface to the model bottom then slices these at land surface to identify the outcropping High Conductivity Pathways (HCPs).
2. 02_Setback distance geology defines the setbacks was previously used to run a quick summary of the permeameter tests.
	-03_Setback flow result plotting summarizes the arrays of HCPs
3. 04_Setback distance manning's and eco is run to create the XS of the river for each setback distance. Calculates example flood depth array.
4. 05_Recharge volume prepares maps of the conductivity at the surface (1-2 m). Calculates an example recharge using the constant flood depth array.
5. 06_Flow passing recharge is the recharge calculation with flow routing between segments based on the flood types from Whipple et al. 2017 and saves information on streamflow, recharge, and depth. This is translated to the script multiprocess_recharge which is used to run the 100 realizations in parallel.
	- multiprocess_recharge is the parallel script to save the key data
6. 06b_recharge_plotting is the current iteration of recharge result processing which summarizes total recharge, looks at streamflow depletion and maps depth
7. 10_method_figures creates a map of the site, example cross-section with flood depths and plots 3D geologic model figures


# Previous Summary of Setback Distance Analysis 

1. Permeameter for velocity rectangular
First step in process was to create permeameters in modflow to force recharge through the domain to identify ideally the high flow cells, although Graham mentioned this may also make some low conductivity facies appear as high flow.
	- Parallel runs of realizations: Script to iterate over all the permeameters, but there is also a python script writtent to be used by the multiprocessing package to use parallel processing. 

2. Setback distance geology: step in process is to read the output from all 100 permeameter runs identify high flow cells.

3. Setback flow result plotting: step is to plot the output from the setback distance geology results.

4. Setback distance manning and eco: step is to create series of cross-sections to estimate the flood depths for the Cosumnes and to iterate across all of the setback distances.

5. Recharge volume: step was to test approximately how much water would be recharged statically by a given flood event, assuming recharge doesn't affect flood. In this step I realized the hydraulic conductivity I used needed a vertical anisotropy or else the volume of water recharged was unreasonable.

6. Muskingum_recharge: step is to combine the setback distance flood depth estimates and the recharge estimates to iterate so that recharge reduces the floow flow in the channel.

7. Connectivity: step was to return to the connectivity analysis because I noticed that some of the 'high flow' cells were in the fine facies group.
