







# +
gel_dir = gwfm_dir+'/UPW_data'
if 'GHB_UZF_WEL_scaling.csv' in os.listdir(model_ws):
    scaling_factors = pd.read_csv(model_ws+'/GHB_UZF_WEL_scaling.csv',delimiter = ',')
else:
    scaling_factors = pd.read_csv(gel_dir+'/GHB_UZF_WEL_scaling.csv',delimiter = ',')
    scaling_factors.to_csv(model_ws+'/GHB_UZF_WEL_scaling.csv')

scaling_factors_all = scaling_factors.copy()
# convert from K (m/s) to K (m/d)
scaling_factors_all.loc[scaling_factors_all.GroupName.isin(['GHB']),'StartValue']*=86400


scaling_factors = scaling_factors.set_index('ParamName')['StartValue']

# +
# def get_layer_from_elev(elev, botm_slice, nlay):
#     """  Return uppermost model layer (0-based) occupied at least partly by some elevation data
#     Parameters
#     ----------
#     elev: 1D array (n) with elevations matching model elevation units
#     botm: 2D array (nlay, n) with layer elevations of model using same x,y locations at elev1D
#     """
#     elev_lay = np.zeros(len(elev))
#     for k in np.arange(0,nlay-1):
#         for j in np.arange(0,len(elev)):
#             if botm_slice[k,j] > elev[j]:
#                 elev_lay[j] = k + 1
#     return(elev_lay.astype(int))


# -

# # SFR

# +
# replaced by version from Oneto Denier model
XSg = XSg_in.copy()

if scenario != 'no_reconnection':
    ## upstream ##
    # od_breach is the sensor location where the breach was made in the levees for flow to leave the river
#     od_breach = fp_grid_xs[fp_grid_xs['Logger Location']=='OD_Excavation'].copy()
    od_breach = OD_locs[OD_locs['Logger Location']=='OD_Excavation'].sjoin_nearest(XSg.reset_index(), how='inner')
#     od_breach.z_min += od_breach.slope*20
#     od_breach.reach -=0.2
    # second reach to assist diversion into floodplain
    # need to adjust elevation so transfer segment from floodplain diversion to stream is positive slope
    od_return = od_breach.copy()
    od_return['Logger Location'] = 'OD_Exc_Return'
    od_return['iseg'] += 1 # adjust iseg to set sorting order
#     od_breach.reach -=0.2

#     od_return.z_min -= od_return.slope*10
    # all segments starting at the floodplain must have 2 added to accomodate the diversion and return flow
    # includes the breach and return
    XSg.loc[XSg.iseg>=od_breach.iseg.values[0], 'iseg'] +=2
    
    ## downstream ##
    # the swale is where flow returns to the channel
    od_swale = OD_locs[OD_locs['Logger Location']=='SwaleBreach_1'].sjoin_nearest(XSg.reset_index(), how='inner')
    # joining to sfr reaches places the swale return in a more realistic row,column
    od_swale = OD_locs[OD_locs['Logger Location']=='SwaleBreach_1'].sjoin_nearest(grid_sfr_in.reset_index(), how='inner')
    #     all segments after the floodplain return must have 1 added to accomodate the return flow
    # XSg.loc[XSg.iseg>=od_swale.iseg.values[0], 'iseg'] += 1
    
    # add reaches for diversion
    XSg = pd.concat((XSg.reset_index(), od_breach, od_return, od_swale)) 
    XSg = XSg.sort_values('iseg')
else:
    XSg = XSg_in.reset_index().copy()

# +
# may not be needed now
# if scenario != 'no_reconnection':
#     # update grid_sfr to account for added segments (4)
#     sfr_add = XSg.loc[~XSg['Logger Location'].isna(), np.append(grid_sfr_in.columns.values,['Logger Location','Site'])]
#     # specify short lengths for identifying (don't think modflow will allow a zero length)
#     sfr_add.length_m = 1
#     # run first by simply re-writing the reach numbers to be longer
#     grid_sfr = pd.concat((grid_sfr_in,sfr_add))
#     # resort by reach number so they are in the correct spot before renumbering
# #     grid_sfr = grid_sfr.sort_values('reach')
# #     grid_sfr.reach = np.arange(1,len(grid_sfr)+1)

# else:
#     grid_sfr = grid_sfr_in.copy()
# -

# # join iseg to grid_sfr
# xs_sfr = grid_sfr.merge(XSg[['row','column', 'iseg']],how='left')
# if scenario=='reconnection':
#     xs_sfr = grid_sfr.merge(XSg[['row','column','Logger Location','Site', 'iseg']],how='left')


# # RCH



# +
# def load_swb_data(strt_date, end_date, param):
#     nper_tr = (end_date-strt_date).days+1
#     # years and array index 
#     # include final year for specifying index then drop
#     years = pd.date_range(strt_date,end_date+pd.DateOffset(years=1),freq='AS-Oct')
#     yr_ind = (years-strt_date).days
#     years = years[:-1]
#     perc = np.zeros((nper_tr, nrow_p,ncol_p))
#     # need separte hdf5 for each year because total is 300MB
#     for n in np.arange(0,len(years)):
#     #     arr = pc[yr_ind[n]:yr_ind[n+1]]
#         fn = join(uzf_dir, 'basic_soil_budget',param+"_WY"+str(years[n].year+1)+".hdf5")
#         with h5py.File(fn, "r") as f:
#             arr = f['array']['WY'][:]
#             perc[yr_ind[n]:yr_ind[n+1]] = arr
#     return(perc)


# -

# load applied water from soil water budget
# it is the ETc scaled by irrigation efficiencies
# AW = load_swb_data(strt_date, end_date, 'applied_water')
# AW_ss = load_swb_data(ss_strt, strt_date-pd.DateOffset(days=1), 'applied_water')


# # WEL

# +
# old version in array format
# simplified ag well layer with just one layer per well
# ag_row, ag_col = np.where(ET_ag.sum(axis=0)>0)
# ag_well_lay = get_layer_from_elev((dem_data-ag_well_depth_arr*0.9)[ag_row, ag_col], 
#                                   botm[:, ag_row, ag_col], m.dis.nlay)
# ag_well_lay.shape, ag_row.shape
# ag_well_lay = pd.DataFrame(np.transpose((ag_row,ag_col, ag_well_lay)), columns=['row','column','layer'])

# +

# iterate over all row, col and get layers for each well based on "screen" 
# ag_well_lay = np.zeros((1,3))
# for i,j in zip(ag_min_lay.reset_index().row,ag_min_lay.reset_index().column):
#     lays = np.arange(ag_min_lay.layer.loc[i,j],ag_max_lay.layer.loc[i,j]+1)
#     ijk = np.rot90(np.vstack((np.tile(i,len(lays)), np.tile(j,len(lays)),lays)))
#     ag_well_lay = np.vstack((ag_well_lay,ijk))
# # delete filler first row
# ag_well_lay = ag_well_lay[1:]
# ag_well_lay = pd.DataFrame(ag_well_lay.astype(int), columns=['row','column','layer'])

# num_ag_layers = (ag_max_lay - ag_min_lay+1).reset_index()
# # divide ET_ag by the number of layers it will go into
# ET_ag_layered = np.copy(ET_ag)
# ET_ag_layered[:,num_ag_layers.row,num_ag_layers.column] /= num_ag_layers.layer.values

# +
# used by the array AW method, not needed when done by field
# adjustments to allow connection with rows,cols with pumping
# row_col = ag_well_lay.loc[:,['row','column']].rename({'row':'rowi','column':'colj'},axis=1)
# ag_well_lay = ag_well_lay.set_index(['row','column'])
# ag_well_lay['rowi'] = row_col.rowi.values
# ag_well_lay['colj'] = row_col.colj.values

# +
# used by the array AW method, not needed when done by field

# create empty dictionary to fill with stress period data
# wel_ETc_dict = {}
# # end date is not included as a stress period, starting at 1st TR spd (2)
# for t in np.arange(0,nper):
#     wel_i, wel_j = np.where(ET_ag[t, :, :]>0)
#     new_xyz = ag_well_lay.loc[list(zip(wel_i,wel_j))] 
# #     wel_ETc = -ET_ag[t-1,wel_i,wel_j]*delr*delr
# # use new row,cols because there are more layers to use
# #     wel_ETc = -ET_ag_layered[t, new_xyz.rowi, new_xyz.colj]*delr*delr
#     wel_ETc = -ET_ag[t, new_xyz.rowi, new_xyz.colj]*delr*delr
#     # ['layer','row','column', 'flux'] are necessary for WEL package
#     spd_ag = np.stack((new_xyz.layer, new_xyz.rowi, new_xyz.colj,wel_ETc),axis=1)
#     # correct by dropping any rows or cols without pumping as some may be added
#     spd_ag = spd_ag[spd_ag[:,-1]!=0,:]
#     spd_all = np.copy(spd_ag)
#     wel_ETc_dict[t] = spd_all
# -
# # HOB plotting


# +

# voi = hob_gpd[hob_gpd.site_code.isin(hobs_long)].node.unique()
# # voi = 22825
# # def obs_sim_node(voi):
# ny = 3
# nx = int(np.ceil(len(voi)/ny))
# fig,ax=plt.subplots(nx,ny,figsize=(12,3*nx))
# for i,voi_n in enumerate(voi):
#     ax_n = ax[int(i / ny), i % ny] if (nx > 1) else ax[i]
#     hob_n = hob_gpd[hob_gpd.node==voi_n]
#     hob_n.reset_index().plot(x='date',y='obs_val',kind='scatter', ax=ax_n,
#                                                   marker='x', s=40, legend=False)
#     hob_n.reset_index().plot(x='date',y='sim_val',kind='scatter', ax=ax_n, 
#                                                   marker='o', s=40, legend=False)
#     # add time series of simulated data to see true peaks
# #     ts_i = pd.DataFrame(hdobj.get_ts(hob_kij[hob_kij.node==voi_n].kij.values[0]),columns=['spd','sim_val']) 
#     ts_i['dt'] = dt_ref.dt
#     ts_i.plot(x='dt',y='sim_val', ax=ax_n, legend=False)

#     ax_n.set_xlabel('')
#     ax_n.set_ylabel('')
#     S_n = format(gel.ss.array[hob_n.layer.iloc[0],hob_n.row.iloc[0],hob_n.column.iloc[0]],'.1e')
#     K_n = format(gel.hk.array[hob_n.layer.iloc[0],hob_n.row.iloc[0],hob_n.column.iloc[0]], '.1e')

#     ax_n.set_title(str(voi_n)+' K'+K_n+' S'+S_n+'\n'+hob_n.site_code.iloc[0])
# ax_n.legend(['Observed','Simulated'])
# ax[0,0].legend(['Observed','Simulated'])

# # fig.text(-0.03, 0.2, 'Head (m)',rotation='vertical',size=26)
# # fig.text(0.35, -0.05, 'Date',size=26)
# fig.tight_layout()

