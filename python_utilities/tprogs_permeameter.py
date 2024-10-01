
# Module with functions to prepare and post-process permeameter runs of TPROGs data
import numpy as np 
import pandas as pd  
import flopy  
import rasterio

from shapely.ops import Polygon
from rasterio.features import shapes, rasterize


def elev_to_tprogs_layers(elev, tprogs_top_elev, tprogs_bot_elev, num_lays):
    """
    function to get the tprogs layers based on the given elevation
    Example
    layer 0 is 80 meters, layer 1 is 79.5 meters, layer -1 is -80 meters
    """
    lay_thick = (tprogs_top_elev - tprogs_bot_elev)/num_lays
    elev_round = np.round((elev) * (1/lay_thick)) / (1/lay_thick) # dem rounded to the layer thickness
    elev_round[elev_round >= tprogs_top_elev] = tprogs_top_elev# any elevation above the top is set to the top
    # subtract the calculated row from top elev divided by layer thickness to get to index 0 at top and index 320 and bottom
    elev_indices = tprogs_top_elev/lay_thick - elev_round*(1/lay_thick) 
    return(elev_indices.astype(int))

def highflow_at_groundsurface(run_ws, flow_percentile):
    ''' take Cell by Cell budget file and finds high flow cells by percentile
    then finds those that outcrop at ground surface'''
    cbb = flopy.utils.CellBudgetFile(run_ws+'/MF.cbc')
    # load velocity in z direction
    extcbb = flopy.utils.postprocessing.get_extended_budget(cbb)
    (qx, qy, qz) = flopy.utils.postprocessing.get_specific_discharge(vectors = extcbb, model=m)
    # convert flow to positive as it is all moving in the downward, -z direction
    # q = qz * -1 # not a good indicator at all
    # much better to use magntiude of velocity vector
    q = np.sqrt(qx**2 + qy**2 + qz**2)
    # split cells into low and high conductivity, based on chosen flow percentile
    q_lay = np.zeros((320, 100,230))
    q_lay[q >= np.percentile(q,flow_percentile)] = 1

    # get high conductivity at ground surface
    q_plt = np.zeros((100,230))
    q_plt[rows,cols] = q_lay[tprogs_lay[rows,cols],rows,cols] 
    return(q_plt)

def high_flow_arr(flow_percentile, str_setbacks):
    tic = time.time()
    # map the high flow cells for all 100 realizations
    hf_tot = np.zeros((100,nrow,ncol))

    for r in np.arange(0,100):
        print('Realization', r, ' time since start ',(time.time()-tic)/60)
        folder = '/realization'+ str(r).zfill(3)+'/'
        run_ws = model_ws+folder
        
        qz_lay = highflow_at_groundsurface(run_ws, flow_percentile)
        # calculate total cells in each setback
        hf_tot[r,:,:] = qz_lay

    # save counted high flow cells to a tsv
    hf_tot_out = np.reshape(hf_tot, (100*nrow,ncol))
    np.savetxt(data_dir+'surface_highflow_by_realization_'+str(flow_percentile)+'.tsv', hf_tot_out,delimiter = '\t')

    toc = time.time()
    print('Total time was', (toc-tic)/60, 'minutes')

def overlay_hf_setback(qz_lay, str_setbacks):
    qz_lay_setback = np.zeros((100,230)) # overlap high flow cells with setback distance
    qz_lay_setback[str_setbacks.astype('bool')] = qz_lay[str_setbacks.astype('bool')]
    #             hf = qz_lay[str_setbacks[n,:,:].astype('bool')]
    return(qz_lay_setback) # hf.sum()




def calc_area_stats(r, s, l, qz_lay, stat_cols):
    ''' for a given realizaiton and setback, find unique cell groups and take summary stats'''
    # The function shapes from rasterio requires uint8 format
    qz_lay_uint = qz_lay.astype(rasterio.uint8)
    # 'Values of False or 0 will be excluded from feature generation'
    out = shapes(qz_lay_uint, mask=qz_lay.astype(bool), connectivity = 8)
    alldata = list(out)
    cell_stats = pd.DataFrame(np.zeros((1,len(stat_cols))), columns=stat_cols)
    num_cells = np.zeros((len(alldata)))
    # iterate over all high flow cell groups
    for i in np.arange(0,len(alldata)):
        # coordinates are in terms of row and column number
        grp_coords = alldata[i][0].get('coordinates')[0]
        # the polygon area corresponds to the number of cells included
        grp_poly = Polygon(grp_coords)
        # grp_poly = MultiPoint(temp) # method to check corners
        num_cells[i] = grp_poly.area
    if len(num_cells)!=0: #if a realization has no cells then leave as zeros
        # calculate statistics for given realization, setback
        cell_stats = cell_stats.assign(Num_Grps = len(num_cells), Mean = num_cells.mean(), Median = np.median(num_cells), 
                                       Min = num_cells.min(),Max = num_cells.max(), Variance = num_cells.var(),
                                       Realization = r, Setback = s)
    if 'Location' in stat_cols:
        cell_stats.Location = l
    return(cell_stats)


def high_flow_count(flow_percentile, str_setbacks, local_str_setbacks):
    tic = time.time()
    # will count total number of cells for each setback distance and for all 100 realizations
    hf_tot = np.zeros((100,len(setbacks)))
    # layer for each local setback
    num_sites = len(np.unique(local_str_setbacks))-1
    hf_tot_local = np.zeros((num_sites, 100,len(setbacks)))
    hf_all = np.zeros((100, 100, 230)) # map high flow for each realization
    
    # dataframe for grouping and area analysis
    stat_cols = ['Num_Grps','Mean','Median','Min','Max','Variance','Realization', 'Setback']
    cell_stats_all = pd.DataFrame(np.zeros((100*len(str_setbacks),len(stat_cols))), columns=stat_cols)
    local_cols = stat_cols+['Location']
    cell_stats_all_local = pd.DataFrame(np.zeros((num_sites*100*len(str_setbacks),len(local_cols))), columns=local_cols)
    
    k=0 # counter 
    kl = 0 # local counter
    for r in np.arange(0,100):
        print('Realization', r, ' time since start ',(time.time()-tic)/60)
        folder = '/realization'+ str(r).zfill(3)+'/'
        run_ws = model_ws+folder
        
        qz_lay = highflow_at_groundsurface(run_ws, flow_percentile)
        hf_all[r,:] = np.copy(qz_lay)
        # complete analysis for regional and local setbacks
        for n in np.arange(0,len(setbacks)):
            # overlay high flow cells with setback distance
            qz_lay_setback = overlay_hf_setback(qz_lay, str_setbacks[n,:,:])
            # calculate total cells in each setback
            hf_tot[r,n] = qz_lay_setback.sum()
            # calculate high flow groups and summary statistics
            cell_stats_all.iloc[k] = calc_area_stats(r,n, 0, qz_lay_setback, stat_cols)
            # iterate over local setbacks
            for l in np.unique(local_str_setbacks)[1:].astype(int):
                arr = np.zeros(local_str_setbacks[n,:,:].shape)
                arr[local_str_setbacks[n,::]==l] = 1
                qz_lay_setback_local = overlay_hf_setback(qz_lay, arr)
                hf_tot_local[l-1,r,n] = qz_lay_setback_local.sum()
                cell_stats_all_local.iloc[kl] = calc_area_stats(r,n,l, qz_lay_setback_local, local_cols)
                kl+=1
            k +=1
    hf_tot_df = pd.DataFrame(hf_tot, columns = setbacks)
    hf_tot_local = np.reshape(hf_tot_local, (num_sites*100,len(setbacks)))
    hf_tot_local_df = pd.DataFrame(hf_tot_local, columns = setbacks)
    hf_all_out = np.reshape(hf_all, (100*100, 230))
    np.savetxt(data_dir+'surface_highflow_by_realization_'+str(flow_percentile)+'.tsv', hf_all_out, delimiter = '\t')
    
    # save counted high flow cells to a csv
    hf_tot_df.to_csv(data_dir+'surface_highflow_by_distance_regional_'+str(flow_percentile)+'.csv', index=False)
    hf_tot_local_df.to_csv(data_dir+'surface_highflow_by_distance_local_'+str(flow_percentile)+'.csv', index=False)
    # save grouping analysis and area statistics
    cell_stats_all.to_csv(data_dir+'surface_highflow_cells_statistics_regional'+str(flow_percentile)+'.csv', index=False)
    cell_stats_all_local.to_csv(data_dir+'surface_highflow_cells_statistics_local'+str(flow_percentile)+'.csv', index=False)

    toc = time.time()
    print('Total time was', (toc-tic)/60, 'minutes')

