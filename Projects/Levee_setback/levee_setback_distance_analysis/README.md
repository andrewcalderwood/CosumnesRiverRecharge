# Summary of Setback Distance Analysis 

## Permeameter for velocity rectangular
First step in process was to create permeameters in modflow to force recharge through the domain to identify ideally the high flow cells, although Graham mentioned this may also make some low conductivity facies appear as high flow.

### Parallel runs of realizations 
Script to iterate over all the permeameters, but there is also a python script writtent to be used by the multiprocessing package to use parallel processing. 

## Setback distance geology
Second step in process is to read the output from all 100 permeameter runs identify high flow cells.

## Setback flow result plotting
Third step is to plot the output from the setback distance geology results.

# Setback distance manning and eco
Fourth step is to create series of cross-sections to estimate the flood depths for the Cosumnes and to iterate across all of the setback distances.

# Recharge volume
Fifth step was to test approximately how much water would be recharged statically by a given flood event, assuming recharge doesn't affect flood. In this step I realized the hydraulic conductivity I used needed a vertical anisotropy or else the volume of water recharged was unreasonable.

# Muskingum_recharge
Sixth step is to combine the setback distance flood depth estimates and the recharge estimates to iterate so that recharge reduces the floow flow in the channel.

# Connectivity
Seventh step was to return to the connectivity analysis because I noticed that some of the 'high flow' cells were in the fine facies group.