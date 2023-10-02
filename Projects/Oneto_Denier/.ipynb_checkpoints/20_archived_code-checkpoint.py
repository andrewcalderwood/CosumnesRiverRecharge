# %% [markdown]
# # Old WEL package
#

# %% WEL package

## from when a gridded format was used based ETc alone

## only needed when in grid format
# remove evapotranspiration in stream cells
# ag_local[:, sfr_rows, sfr_cols] = 0

# # ET_ag = np.copy(ag_local)
# ET_ag = np.copy(AW_local)

# if ss_bool == True:
#     # ET_ag_SS = np.reshape(ss_ag_local.mean(axis=0),(1, nrow,ncol))
#     ET_ag_SS = np.copy(AW_ss_local)
#     ET_ag = np.concatenate((ET_ag_SS, ET_ag), axis=0)
# print(ET_ag[0].sum().round(1))
# # remove pumping where it is considered GDE, reduced SS by almost 1/2 (5.7 to 3.6)
# # ET_ag[:, (GDE_cell.row-1).astype(int), (GDE_cell.column-1).astype(int)] = 0
# # print(ET_ag[0].sum().round(1))
# # remove pumping where it is known restoration (floodplain), doesn't remove much (3.65 to 3.61)
# ET_ag[:, lakarr[0] >0] = 0
# print(ET_ag[0].sum().round(1))
# # no pumping where land is fallows
# ET_ag[:, fallow.row-1, fallow.column-1] = 0
# print(ET_ag[0].sum().round(1))

# removes 80% of wells
# ag_screen_botm
# remove wells that have more than 20% below the model bottom
# drop_wells = pd.DataFrame(np.rot90(np.where((dem_data-ag_well_depth_arr)*0.8 < botm[-1])), columns=['row','column'])
# drop_wells = pd.DataFrame(np.rot90(np.where((dem_data-ag_well_depth_arr)*0.8 < -37)), columns=['row','column'])

# ag_well_elev = dem_data - ag_well_depth_arr
# ag_screen_botm = np.where((ag_well_elev<botm)&(ag_well_elev> botm[-1]))
# ag_screen_botm = np.rot90(ag_screen_botm)
# ag_screen_botm = pd.DataFrame(ag_screen_botm, columns=['layer','row','column'])
# ag_max_lay = ag_screen_botm.groupby(['row','column']).max()
# # any wells below most bottom go in bottom layer
# ag_max_lay.layer[ag_max_lay.layer == nlay] = nlay-1

# # assume 10% of well is screened? Pauloo? tprogs lay thickness is 4m, so 12ft, not quite enough for typical well?
# # if we go two layers we have 8 m which is near the average expected well screen
# ag_screen_top = np.where((ag_well_elev*0.8 <botm)&(ag_well_elev*0.8>botm[-1]))
# ag_screen_top = np.rot90(ag_screen_top)
# ag_screen_top = pd.DataFrame(ag_screen_top, columns=['layer','row','column'])
# ag_min_lay = ag_screen_top.groupby(['row','column']).max()
# ag_min_lay.layer[ag_min_lay.layer == nlay] = nlay-1


# ag_well_lay.layer.median() 
# mean layer is 10, median is 11
# the issue with pumping could be so much in the deeper aquifer it causes issues
# ag_min_lay.shape, ag_max_lay.shape


# # if the layer max is missing set as model bottom
# ag_lays = ag_min_lay.join(ag_max_lay, rsuffix='_mx',lsuffix='_mn')
# ag_lays.loc[ag_lays.layer_mx.isna(),'layer_mx'] = m.dis.nlay-1
# ag_lays.layer_mx = ag_lays.layer_mx.astype(int)

# # iterate over all row, col and get layers for each well based on "screen" 
# ag_well_lay = np.zeros((1,3))
# for i,j in zip(ag_min_lay.reset_index().row,ag_min_lay.reset_index().column):
#     lays = np.arange(ag_lays.layer_mn.loc[i,j], ag_lays.layer_mx.loc[i,j]+1)
#     ijk = np.rot90(np.vstack((np.tile(i,len(lays)), np.tile(j,len(lays)),lays)))
#     ag_well_lay = np.vstack((ag_well_lay,ijk))
# # delete filler first row
# ag_well_lay = ag_well_lay[1:]
# ag_well_lay = pd.DataFrame(ag_well_lay.astype(int), columns=['row','column','layer'])

# num_ag_layers = (ag_max_lay - ag_min_lay+1).reset_index()
# # divide ET_ag by the number of layers it will go into
# ET_ag_layered = np.zeros(ET_ag.shape)
# # ET_ag_layered = np.copy(ET_ag)
# ET_ag_layered[:,num_ag_layers.row,num_ag_layers.column] = ET_ag[:,num_ag_layers.row,num_ag_layers.column]/num_ag_layers.layer.values
# # adjustments to allow connection with rows,cols with pumping
# row_col = ag_well_lay.loc[:,['row','column']].rename({'row':'rowi','column':'colj'},axis=1)
# ag_well_lay = ag_well_lay.set_index(['row','column'])
# ag_well_lay['rowi'] = row_col.rowi.values
# ag_well_lay['colj'] = row_col.colj.values

# layer for ETc ag well pumping
# ETc_lay = 1
# create empty dictionary to fill with stress period data
# wel_ETc_dict = {}
# # end date is not included as a stress period, starting at 1st TR spd (2)
# for t in np.arange(0,nper):
#     wel_i, wel_j = np.where(ET_ag_layered[t,:,:]>0)
#     new_xyz = ag_well_lay.loc[list(zip(wel_i,wel_j))] 
# # use new row,cols because there are more layers to use
#     wel_ETc = -ET_ag_layered[t,new_xyz.rowi,new_xyz.colj]*delr*delr
#     # ['layer','row','column', 'flux'] are necessary for WEL package
#     spd_ag = np.stack((new_xyz.layer, new_xyz.rowi, new_xyz.colj,wel_ETc),axis=1)
#     # correct by dropping any rows or cols without pumping as some may be added
#     spd_ag = spd_ag[spd_ag[:,-1]!=0,:]
# #     dom_loc['flux'] = - dom_use.loc[dates[t],'flux_m3d']
# #     wells_dom = dom_loc[['layer','row','column','flux']].values
# #     spd_noag = np.vstack((wells_dom))
#     # join pumping from ag with point pumping from domstic/supply wells that are constant
# #     spd_all = np.vstack((spd_ag,spd_noag)) 
#     spd_all = np.copy(spd_ag)
#     wel_ETc_dict[t] = spd_all

# %% [markdown]
# # Post-processing

# %%
# hdobj = flopy.utils.HeadFile(model_ws+'/MF.hds')
# # extract time series of heads for each desired location
# mw_hds = hdobj.get_ts(list(zip(rm_grid['lay'], hob_row, hob_col)))
# mw_hds = pd.DataFrame(mw_hds, columns=['time']+rm_grid.Sensor.tolist())
# # convert to hourly to maintain more precision in DT
# mw_hds['dt'] = strt_date+(mw_hds.time.values*24 ).astype('timedelta64[h]')
# mw_gwl = mw_hds.drop(columns=['time'])
# # long format for id join with observed dat
# mw_long = mw_gwl.melt(id_vars='dt', var_name='Well',value_name='sim')
# mw_long = mw_long[mw_long.sim != -1e30]

# %%
# mw_chk = mw_long.join(gwl_long.set_index(['Well','dt']), on=['Well','dt'], how='inner')
# mw_chk = mw_chk.melt(id_vars=['dt', 'Well'],value_vars=['sim','obs'], value_name='gwe', var_name='type')

# %%
# sns.relplot(mw_chk,x='dt',y='gwe',col='Well', hue='type', col_wrap=4)

# %%

# %%
