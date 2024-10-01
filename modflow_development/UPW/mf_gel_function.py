
# %% [markdown]
# ## LPF/UPW package

# %%
hk = np.zeros(botm.shape)
vka = np.zeros(botm.shape)
sy = np.zeros(botm.shape)
ss = np.zeros(botm.shape)

# %%
top = np.copy(botm[0,:,:]) # bottom of levelling layer
bot1 = np.copy(botm[-3,:,:]) # top of laguna

# I need to verify if a flattening layer is needed (e.g., variable thickness to maintain TPROGs connectivity)
# pull out the TPROGS data for the corresponding depths
K_c = tc.get_tprogs_for_elev(K, top, bot1, tprogs_info)
Ss_c = tc.get_tprogs_for_elev(Ss, top, bot1, tprogs_info)
Sy_c = tc.get_tprogs_for_elev(Sy, top, bot1, tprogs_info)
n_c = tc.get_tprogs_for_elev(porosity, top, bot1, tprogs_info)

tprog_strt_lay = num_leveling_layers+drain_layer
# upscale as preset
for kt, k in enumerate(np.arange(tprog_strt_lay,nlay_tprogs+tprog_strt_lay)):
    hk[k,:] = np.mean(K_c[upscale*kt:upscale*(kt+1)], axis=0)
    vka[k,:] = hmean(K_c[upscale*kt:upscale*(kt+1)], axis=0)
    ss[k,:] = np.mean(Ss_c[upscale*kt:upscale*(kt+1)], axis=0)
    sy[k,:] = np.mean(Sy_c[upscale*kt:upscale*(kt+1)], axis=0)
#     por[k,:] = np.mean(n_c[upscale*kt:upscale*(kt+1)], axis=0)

# %%
top = m.dis.top.array
bot1 = m.dis.botm.array[drain_layer,:,:]
# set parameters based on upscaled unsaturated zone
for k in np.arange(0, tprog_strt_lay):
    hk[k,:,:] = np.mean(tc.get_tprogs_for_elev(K, top, bot1,tprogs_info),axis=0)
    vka[k,:,:] = hmean(tc.get_tprogs_for_elev(K, top, bot1,tprogs_info),axis=0)
    sy[k,:,:] = np.mean(tc.get_tprogs_for_elev(Sy, top, bot1,tprogs_info),axis=0)
    ss[k,:,:] = np.mean(tc.get_tprogs_for_elev(Ss, top, bot1,tprogs_info),axis=0)



# %%
rows,cols = grid_p.row.values-1, grid_p.column.values-1

# %%
# check proportions of hydrofacies in TPROGs realization
tprogs_vals = np.arange(1,5)
tprogs_hist = np.histogram(masked_tprogs, np.append([0],tprogs_vals+0.1))[0]
tprogs_hist = tprogs_hist/np.sum(tprogs_hist)

tprogs_quants = 1 - np.append([0], np.cumsum(tprogs_hist)/np.sum(tprogs_hist))
vka_quants = pd.DataFrame(tprogs_quants[1:], columns=['quant'], index=tprogs_vals)
# dataframe summarizing dominant facies based on quantiles
vka_quants['vka_min'] = np.quantile(vka, tprogs_quants[1:])
vka_quants['vka_max'] = np.quantile(vka, tprogs_quants[:-1])
vka_quants['facies'] = params.loc[tprogs_vals].Lithology.values
# scale vertical conductivity with a vertical anisotropy factor based
# on quantiles in the upscaled tprogs data
for p in tprogs_vals:
    vka[(vka<vka_quants.loc[p,'vka_max'])&(vka>vka_quants.loc[p,'vka_min'])] /= params.vani[p]


# %%
vka_quants.to_csv(join(model_ws, 'vka_quants.csv'))

# %% [markdown]
# The tuff breccia is very dense, hard and low water yielding. It is supposedly responsible for the many "haystack" hills in the eastern part of the county
#
# DWR report has a few final well pumping rates, drawdowns and specific capacities but limited.
#
# Fleckenstein et al. 2006 found the Mehrten had  
# Kh = 1 to 1.8 x10^-5 m/s  
# Kv = 1 to 1.8 x10^-7 m/s  
# vani ~= 100  
# Sy = 0.15 to 0.2  
# Ss = 1e-4 to 1e-3 m^-1  

# %%
# set values for second to bottom layer, Laguna formation
hk[-2,:,:] = params.loc[5,'K_m_d']
vka[-2,:,:] = params.loc[5,'K_m_d']/params.loc[5,'vani'] 
sy[-2,:,:] = params.loc[5,'Sy']
ss[-2,:,:] = params.loc[5,'Ss']

# set values for bottom layer, Mehrten formation
hk[-1,:,:] = params.loc[6,'K_m_d']
vka[-1,:,:] = params.loc[6,'K_m_d']/params.loc[6,'vani'] 
sy[-1,:,:] = params.loc[6,'Sy']
ss[-1,:,:] = params.loc[6,'Ss']

# %%
# rather than use a variable deep geology array which is complicated to determine local effects
# use the mean column for each layer to define a block of Low K to correct gradient in the foothill
adj_lowK = pd.DataFrame(np.transpose(np.where(deep_geology>0)), columns=['k','i','j'])
# the mean didn't quite extend far enough or wasn't low enough K
# adj_lowK = adj_lowK.groupby('k').mean()['j'].astype(int)
# trying near minimum to extend further, manually adjusted to 0.15 by comparing to dem_data>56
adj_lowK = adj_lowK.groupby('k').quantile(0.15)['j'].astype(int)
adj_lowK_arr = np.zeros((nlay,nrow,ncol))
for k in adj_lowK.index:
    adj_lowK_arr[k, :, adj_lowK.loc[k]:] = 1
# convert to boolean
adj_lowK_arr = adj_lowK_arr.astype(bool)
# don't want to adjust deepest two layers?
# this doesn't make as much sense geologically
# adj_lowK_arr[-1] = 0
# adj_lowK_arr[-2:] = 0
# this is causing potentially high water levels in the foothills
# the deep_geology array shows where the mehrten formation comes out of the surface
# param_d = 7 # distinct parameter for foothills
param_d = 6 # reuse Mehrten parameter (simpler)
hk[adj_lowK_arr] = params.loc[param_d,'K_m_d']
vka[adj_lowK_arr] = params.loc[param_d,'K_m_d']/params.loc[param_d,'vani'] # quick test of if only 10 vani (fractured rock has lower vani)
sy[adj_lowK_arr] = params.loc[param_d,'Sy']
ss[adj_lowK_arr] = params.loc[param_d,'Ss']

# if drain layer is active then use higher K Laguna for top layer
if drain_layer==1:
    hk[0,adj_lowK_arr[0]] = params.loc[5,'K_m_d']
    vka[0,adj_lowK_arr[0]] = params.loc[5,'K_m_d']/params.loc[5,'vani']