''' dis functions module
First iteraton: 2024-8-2
'''

import numpy as np

def set_layer_botm(botm,dem_data, 
                    tprog_strt, tprog_thick, nlay_tprogs,
                   laguna_bot, mehrten_bot,
                   num_leveling_layers, drain_layer):
    botm[num_leveling_layers-1,:,:] = tprog_strt
    # for i in np.arange(num_leveling_layers-1,0,-1):
    #     botm[i-1,:,:] = botm[i,:,:] - leveling_layer_thickness[i]
     
    # if a drain layer is active then move leveling layer down
    if drain_layer ==1:
        botm[1] = tprog_strt
        botm[0] = dem_data - 1
        # drain layer must be above 2nd layer
        botm[0] = np.where(botm[0] <= botm[1] + 1 , botm[1] + 1, botm[0])
       
    # Create TPROGS layers from top down
    tprog_strt_lay = num_leveling_layers+drain_layer
    for i in np.arange(tprog_strt_lay, tprog_strt_lay + nlay_tprogs):
        botm[i,:,:] = botm[i-1,:,:] -tprog_thick
        
    # Thickness to give to bottom layers below the TPROGS layers just to provide adequate spacing,
    # this will be corrected by changing the geology in the layers above to account for what is actually in
    # the Mehrten and what is in the Laguna formations, thickness of 5 also prevents any messy overlap
    thickness_to_skip =10
    botm[-2] = laguna_bot
    # need to adapat this for the case when the bottom of a layer doesn't impact the laguna
    # # Find where top boundary of Mehrten Formation rises within 10 meters of the top layer (10m for sufficient layer thickness)
    bot3ind = np.min(np.where(botm[-2,:,:]>botm[-3,:,:]- thickness_to_skip)[1])
    
    # # Where the top boundary of Mehrten was within 10 meters of the top layer 
    # # set it equal to top layer elevation minus 10 for sufficient layer thickness
    botm[-2,:,bot3ind:] = botm[-3,0,bot3ind]- thickness_to_skip
    botm[-1] = mehrten_bot
    # # Repeat steps above for bottom of Mehrten formation with the top of the Mehrten formation
    bot3ind = np.min(np.where(botm[-1,0,:]>botm[-2,0,:]- thickness_to_skip))
    botm[-1,:,bot3ind:] = botm[-2,0,bot3ind]-thickness_to_skip
    return(botm)

