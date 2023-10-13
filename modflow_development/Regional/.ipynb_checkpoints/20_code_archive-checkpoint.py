







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


