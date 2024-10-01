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

# %% [markdown]
# ## Load WB

# %%
# switch to module function 10/12/2023
# def clean_wb(model_ws, dt_ref, name='flow_budget.txt'):
#     # load summary water budget
#     wb = pd.read_csv(join(model_ws, name), delimiter=r'\s+')
#     # wb = pd.read_csv(loadpth+'/oneto_denier_upscale8x_2014_2018'+'/flow_budget.txt', delimiter=r'\s+')

#     wb['kstpkper'] = list(zip(wb.STP-1,wb.PER-1))
#     wb = wb.merge(dt_ref, on='kstpkper')
#     wb = wb.set_index('dt')
#     # calculate change in storage
#     wb['dSTORAGE'] = wb.STORAGE_OUT - wb.STORAGE_IN
#     # calculate the cumulative storage change
#     wb['dSTORAGE_sum'] = wb.dSTORAGE.cumsum()
#     # calculate net groundwater flow
#     wb['GHB_NET'] = wb.GHB_IN - wb.GHB_OUT
#     # identify relevant columns
#     wb_cols = wb.columns[wb.columns.str.contains('_IN|_OUT')]
#     wb_cols = wb_cols[~wb_cols.str.contains('STORAGE')]
#     wb_out_cols= wb_cols[wb_cols.str.contains('_OUT')]
#     wb_in_cols = wb_cols[wb_cols.str.contains('_IN')]
#     # only include columns with values used
#     wb_out_cols = wb_out_cols[np.sum(wb[wb_out_cols]>0, axis=0).astype(bool)]
#     wb_in_cols = wb_in_cols[np.sum(wb[wb_in_cols]>0, axis=0).astype(bool)]
#     return(wb, wb_out_cols, wb_in_cols)

# %% [markdown]
# ## t-test

# %%

# switched to multiple functions
# def run_stats(wb, wb0, term, season=None, freq='monthly', plot=False):
#     if season == 'Wet':
#         wet_months = [11,12,1,2,3,4]
#         wb = wb[wb.index.month.isin(wet_months)]
#         wb0 = wb0[wb0.index.month.isin(wet_months)]
#     elif season =='Dry':
#         dry_months = [5,6,7,8,9,10]
#         wb = wb[wb.index.month.isin(dry_months)]
#         wb0 = wb0[wb0.index.month.isin(dry_months)]
#     elif season=='Fall':
#         fall_months=[9,10,11]
#         wb = wb[wb.index.month.isin(fall_months)]
#         wb0 = wb0[wb0.index.month.isin(fall_months)]
#     if freq=='annual':
#         a = wb0.resample('AS-Oct').sum()[term].values
#         b = wb.resample('AS-Oct').sum()[term].values
#     elif freq=='monthly':
#         a = wb0.resample('MS').sum()[term].values
#         b = wb.resample('MS').sum()[term].values
#     elif freq=='daily':
#         a = wb0.resample('D').sum()[term].values
#         b = wb.resample('D').sum()[term].values
        
#     t_out = ttest_rel(a, b)

#     t_df = pd.DataFrame([t_out.statistic, t_out.pvalue]).transpose()
#     t_df.columns=['statistic','pvalue']
#     t_df['term'] = term
#     t_df['freq'] = freq
#     t_df['season'] = season
#     t_df['mean_a'] = np.mean(a)
#     t_df['mean_b'] = np.mean(b)
#     t_df['perc_diff_in_means'] = 100*(np.mean(a)-np.mean(b))/np.abs((np.mean(a)+np.mean(b))/2)

#     # rounding to clean up output
#     t_df.statistic = t_df.statistic.round(3)
#     t_df.pvalue = t_df.pvalue.round(4)
#     t_df.perc_diff_in_means = t_df.perc_diff_in_means.round(2)

#     # if pvalue is insignificant then don't include
#     t_df.loc[t_df.pvalue>=0.05,'perc_diff_in_means'] = '-'
    
#     if plot:
#         slope, intercept, r_value, p_value, std_err = linregress(a, b)
#         print('T-test statistic: %.2f' %t_out.statistic, 'and pvalue: %.4f' %t_out.pvalue)

#         plt.scatter(a, b)
#         x_range = np.array([[np.min((a,b))], [np.max((a,b))]])
#         plt.plot(x_range, slope*x_range + intercept, color='black', linewidth=1)
#         plt.annotate('y = '+str(np.round(slope,3))+'x + '+ str(np.round(intercept,2)), (0.1,0.8), xycoords='axes fraction')
#         plt.title(term)
#         plt.ylabel('No Reconnection')
#         plt.xlabel('Baseline')
#     return(t_df)

# %% [markdown]
# ## Load SFR

# %%
# grid_sfr = pd.DataFrame().from_records(m.sfr.reach_data).rename(columns={'i':'row','j':'column'})
# grid_sfr[['row','column']] += 1 # convert to 1 based to match with SFR output
# pd_sfr = grid_sfr.set_index(['iseg','ireach'])[['rchlen','strtop', 'facies']]
# pd_sfr['Total distance (m)'] = pd_sfr['rchlen'].cumsum()

# switch to module function 10/12/2023
# def clean_sfr_df(model_ws, pd_sfr=None):
#     sfrout = flopy.utils.SfrFile(join(model_ws, m.name+'.sfr.out'))
#     sfrdf = sfrout.get_dataframe()
#     sfrdf = sfrdf.join(dt_ref.set_index('kstpkper'), on='kstpkper').set_index('dt')
#     # convert from sub-daily to daily using mean, lose kstpkper
#     sfrdf = sfrdf.groupby('segment').resample('D').mean(numeric_only=True)
#     sfrdf = sfrdf.reset_index('segment', drop=True)
#     sfrdf[['row','column']]-=1 # convert to python
#     sfrdf['month'] = sfrdf.index.month
#     sfrdf['WY'] = sfrdf.index.year
#     sfrdf.loc[sfrdf.month>=10, 'WY'] +=1
#     # add column to track days with flow
#     sfrdf['flowing'] = 1
#     sfrdf.loc[sfrdf.Qout <= 0, 'flowing'] = 0
#     if pd_sfr is not None:
#     #     sfrdf = pd_sfr.join(sfrdf.set_index(['row','column']),on=['row','column'],how='inner',lsuffix='_all')
#         sfrdf = sfrdf.join(pd_sfr ,on=['segment','reach'],how='inner',lsuffix='_all')

#     # create different column for stream losing vs gaining seeapge
#     sfrdf['Qrech'] = np.where(sfrdf.Qaquifer>0, sfrdf.Qaquifer,0)
#     sfrdf['Qbase'] = np.where(sfrdf.Qaquifer<0, sfrdf.Qaquifer*-1,0 )
#     # booleans for plotting
#     sfrdf['gaining'] = (sfrdf.gradient == 0)
#     sfrdf['losing'] = (sfrdf.gradient >= 0)
#     sfrdf['connected'] = (sfrdf.gradient < 1)
#     return(sfrdf)



# %%

# %% [markdown]
# # Plotting

# %% [markdown]
# ## Seepage and flow time series

# %%

# %% [markdown]
# ## Correlations

# %%
# yearly
# the data is already summed across the year and should be averaged across segments
# corr_bool = sfr_yr_sum_all.groupby(['realization', 'WY']).mean()
# corr_bool = corr_bool[['flowing','connected','gaining']]
# # seepage data just needs to be summed again 
# corr_seep = sfr_yr_sum_all.groupby(['realization', 'WY']).sum(numeric_only=True)[['Qbase','Qrech']]
# # flow data should be averaged
# corr_flow = sfr_last_all.groupby(['realization','WY']).mean(numeric_only=True)[['Qout']]
# # join together the data for correlations
# corr_all = corr_seep.join(corr_bool).join(corr_flow).reset_index('WY')

# %%
# corr_out = pd.DataFrame()
# for n_wy in wy_vals:
# # corr_all.groupby('WY').apply(lambda x : pearsonr(coarse_ref.num_coarse, x))
#     corr_wy = corr_all[corr_all.WY==n_wy].drop(columns=['WY'])
#     corr_out = pd.concat((corr_out, calc_corr_stats(corr_wy, coarse_ref.num_coarse).assign(WY=n_wy)))

# %%
# fig, ax = plt.subplots(3,2, figsize=(8,5.3), sharex=True, sharey=True, layout='constrained')
# for n, v in enumerate(variables.keys()):
#     ax_n = ax[int(n/2), n%2]
#     corr_plt = corr_out.loc['r'][[v,'type','WY']].pivot_table(index='WY', values=v, columns='type')
#     # correct order of tests which were resorted
#     corr_plt[tests].plot(kind='bar', ax=ax_n, legend=False, rot=0)
#     ax_n.set_title(variables[v])
# ax[0,1].legend(loc='best')

# %% [markdown]
# When aggregated to a yearly level the correlations don't really change between water years.

# %%

# %%
