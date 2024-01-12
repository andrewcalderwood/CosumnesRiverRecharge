# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import os
import sys
from os.path import join, exists, basename, dirname
import glob
from importlib import reload

import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes, mark_inset
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# +
doc_dir = os.getcwd()
while basename(doc_dir) != 'Documents':
    doc_dir = dirname(doc_dir)
# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'
# dir of stream level data for seepage study
proj_dir = join(gwfm_dir, 'Oneto_Denier')
dat_dir = join(proj_dir,'Stream_level_data')
map_dir = join(gwfm_dir,'Mapping')

sfr_dir = join(gwfm_dir,'SFR_data')
hob_dir = join(gwfm_dir,'HOB_data')
uzf_dir = join(gwfm_dir, 'UZF_data')

# -

out_dir = join(proj_dir, 'output')
fig_dir = join(proj_dir, 'figures')

# +
from importlib import reload
def add_path(fxn_dir):
    if fxn_dir not in sys.path:
        sys.path.append(fxn_dir)
        
flopy_dir = doc_dir+'/GitHub/flopy'
if flopy_dir not in sys.path:
    sys.path.insert(0, flopy_dir)

import flopy 
import flopy.utils.binaryfile as bf


add_path(doc_dir+'/GitHub/CosumnesRiverRecharge/tprogs_utilities')
add_path(doc_dir+'/GitHub/CosumnesRiverRecharge/python_utilities')
# -

import map_cln
reload(map_cln)
from map_cln import gdf_bnds, pnt_2_tup, lab_pnt, plt_cln, arr_lab, xy_lab, make_multi_scale, dir_arrow, plt_arrow
from mf_utility import get_dates, get_layer_from_elev, clean_wb
from flopy_utilities import reach_data_gdf
from tprogs_review import get_tprogs_quants

# +
# import mf_utility
# reload(mf_utility)

# +
# fig, ax = plt.subplots()
# m_domain.plot(ax=ax, color='none', edgecolor='black')
# plt_arrow(ax)
# plt_scale(ax)

# plt_arrow(ax, 0.05, 0.15)
# make_multi_scale(ax, 0.1, 0.1, dist = 0.5E3, scales = [4,2,1])
# # make_multi_scale(ax, 0.1, 0.1, dist = 1E3, scales = [3,2,1])

# -

dist_lab = str(np.round(0.5, 1)).replace('.0','')
len(dist_lab)

grid_r = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')
grid_r = grid_r.rename(columns={'node':'node_p', 'row':'row_p','column':'column_p'})


# +
# write modelgrid to get updated row and col numbers specific to the child grid
grid_dir = join(gwfm_dir, 'DIS_data/streambed_seepage/grid')
grid_fn = join(grid_dir, 'inset_oneto_denier','rm_only_grid.shp')

# m.modelgrid.write_shapefile(grid_fn)
grid_p = gpd.read_file(grid_fn)
grid_p.crs = 'epsg:32610'
# -

rot = 52.9

# the spatial join has duplicate row,columns of the local
# grid_match = gpd.sjoin(grid_p, grid_r, how='left')
grid_match = gpd.overlay(grid_p, grid_r)
grid_match = grid_match.loc[grid_match.geometry.area>9999].reset_index() # correct any duplicates due to tiny geometry inconsistencies

grid_p['id'] = 0
m_domain = grid_p.dissolve('id')

# get exterior polyline of model grid
grid_bnd = gpd.GeoDataFrame(pd.DataFrame([0]), geometry = [grid_p.unary_union.exterior], crs=grid_p.crs)
# find cells that construct the model boundary
bnd_cells = gpd.sjoin(grid_p, grid_bnd)
# bnd_cells = bnd_cells[['row','column']] - 1
# bnd_cells['grid_id'] = np.arange(0,len(bnd_cells))
# bnd_rows, bnd_cols = bnd_cells.row.values, bnd_cells.column.values

# +
rivers = gpd.read_file(join(sfr_dir,"Sac_valley_rivers/Sac_valley_rivers.shp"))
rivers = rivers.to_crs('EPSG:32610')

mb_regional = gpd.read_file(join(gwfm_dir,"DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp"))
# mb = gpd.read_file(join(gwfm_dir,"DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp"))
rivers_clip = gpd.clip(rivers, mb_regional)

soam = gpd.read_file(join(map_dir,"so_am_subbasin/so_am_subbasin.shp"))
soam = soam.to_crs('EPSG:32610')
cos = gpd.read_file(join(map_dir,"cos_subbasin/cos_subbasin.shp"))
cos = cos.to_crs('EPSG:32610')

ca = gpd.read_file(join(map_dir,"ca_state_boundary/CA_State_TIGER2016.shp"))
ca = ca.to_crs('EPSG:32610')
# -

cr = rivers_clip[rivers_clip.GNIS_Name=='Cosumnes River']
mr = rivers_clip[rivers_clip.GNIS_Name=='Mokelumne River']

prms_gis = join(gwfm_dir,'PRMS','GIS')
ws_gdf = gpd.read_file(join(prms_gis, "NHD_H_18040013_HU8_Shape\Shape\WBDHU8.shp"))
ws_gdf=ws_gdf.to_crs(m_domain.crs)

# +
df = pd.read_csv(join(hob_dir,'CosumnesRiverPreserve_MW_screened_interval.csv'))
rm_sp = gpd.GeoDataFrame(df, geometry = gpd.points_from_xy(df.Longitude,df.Latitude), crs='epsg:4326')
rm_sp = rm_sp.to_crs('epsg:32610')
rm_sp = rm_sp.rename(columns={'Well ID':'Sensor'})

# prepare output for modelgrid join
rm_t = rm_sp[rm_sp['At Oneto-Denier']=='Yes']

# -

levee_removal = gpd.read_file(join(map_dir, 'levee_removal_OD','levee_removal_OD.shp'))

# +
lak_shp = join(gwfm_dir,'LAK_data/floodplain_delineation')
# ifp = gpd.read_file(join(lak_shp,'inner_floodplain_domain/inner_floodplain_domain.shp' )).to_crs('epsg:32610')
# lfp = gpd.read_file(join(lak_shp,'lower_floodplain_approximate_area/lower_floodplain_approximate_area.shp' )).to_crs('epsg:32610')
lak_extent = gpd.read_file(join(lak_shp,'LCRFR_ModelDom_2017/LCRFR_2DArea_2015.shp' )).to_crs('epsg:32610')

fp_logger = pd.read_csv(join(gwfm_dir,'LAK_data','floodplain_logger_metadata.csv'))
fp_logger = gpd.GeoDataFrame(fp_logger, geometry = gpd.points_from_xy(fp_logger.Easting, fp_logger.Northing), crs='epsg:32610')
# find grid cell it is within
# fp_grid = gpd.sjoin(fp_logger, grid_p, how='left',op='within')
# -

od_breach = fp_logger[fp_logger['Logger Location']=='OD_Excavation']
od_swale = fp_logger[fp_logger['Logger Location']=='SwaleBreach_1']

# +


map_lgd = [
    Patch(facecolor='none', edgecolor='black', linestyle = '--', alpha=1., label='Watershed Extent'),
    Patch(facecolor='none', edgecolor='black',alpha=1., label='Model Extent'),
    Patch(facecolor='red', edgecolor='red',alpha=0.5, label='Levee removal and\nexcavation areas'),
    # Patch(facecolor='tab:blue', alpha=0.5, label='Reconnected Floodplain'),
    Patch(facecolor='none',hatch='///', edgecolor='blue',  alpha=0.7, label='Reconnected Floodplain'),
    Line2D([0], [0],color='tab:blue',label='Cosumnes River'),
    Line2D([0], [0], marker='.', linestyle='', color='black', label='Monitoring Well'),
]
lgd_short = [
    # Patch(facecolor='none', edgecolor='black', linestyle = '--', alpha=1., label='Watershed Extent'),
    # Patch(facecolor='none', edgecolor='black',alpha=1., label='Model Extent'),
    Patch(facecolor='red', edgecolor='red',alpha=0.5, label='Levee removal and\nexcavation areas'),
    # Patch(facecolor='tab:blue', alpha=0.5, label='Reconnected Floodplain'),
    # Line2D([0], [0],color='tab:blue',label='Cosumnes River'),
    Line2D([0], [0], marker='.', linestyle='', color='black', label='Monitoring Well'),
]


# +

def ref_map(axins, fontsize=10, lw=1):
    ca.plot(ax = axins,alpha = 0.2, edgecolor='black')
    axins.annotate(text='California', xy=lab_pnt(ca), 
                xytext = (6,-30), textcoords = 'offset points', fontsize=fontsize, 
                bbox=dict(boxstyle="square,pad=0.3", fc="lightgrey", ec="black", lw=lw), zorder=1)
    # second one is smaller inset
    axins.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)

    ws_gdf.plot(color='none',ax=axins, linestyle='--')
    # ws_gdf.plot(color='none',ax=axins, linestyle='--')
    m_domain.plot(color="none",edgecolor='black',ax=axins)
    return axins



# +
fig, axins =plt.subplots(figsize=(6.5,6.5))

# second one is smaller inset
axins.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)

ca.plot(ax = axins,alpha = 0.2, edgecolor='black')
axins.annotate(text='California', xy=lab_pnt(ca), 
            xytext = (6,6), textcoords = 'offset pixels',
            bbox=dict(boxstyle="square,pad=0.3", fc="lightgrey", ec="black", lw=2))
ws_gdf.plot(color='none',ax=axins, linestyle='--')
m_domain.plot(color="none",edgecolor='black',ax=axins)

# first one is CA map
ax = inset_axes(axins, width="100%", height="100%", bbox_to_anchor=(0.8, -.1, 1.3, 1.3),
                  bbox_transform=axins.transAxes, loc=2)
# ax.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)
def plt_arrow(ax, xoff = 0.7, yoff=0.15):
    x, y, arrow_length = xoff, yoff, 0.1
    ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
                arrowprops=dict(facecolor='black', width=5, headwidth=15),
                ha='center', va='center', fontsize=20, 
                xycoords=ax.transAxes)
def plt_scale(ax):
    fontprops = fm.FontProperties(size=18)
    scalebar = AnchoredSizeBar(ax.transData,
                               2000, '2 km', 'lower right', pad=0.3, sep=2,color='black',
                               frameon=False, size_vertical=1E2, fontproperties=fontprops)
    ax.add_artist(scalebar)
def main_map(ax):
    ## main figure ##
    ax.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)

    plt_bnds = gdf_bnds(m_domain,ax=ax, buf=200)
    # plt_bnds.plot(ax=axins, color='none')
    
    # od_breach.plot(ax=ax)
    # od_swale.plot(ax=ax)
    
    m_domain.plot(color="none", edgecolor='black',ax=ax)
    
    cr.plot(ax=ax, color='tab:blue')
    mr.plot(ax=ax, color='tab:blue')
    
    # lak_extent.plot(ax=ax, color='tab:blue',  alpha=0.5)
    lak_extent.plot(ax=ax, color='none',hatch='///', edgecolor='white',  alpha=0.7)
    # rm_t.plot(legend=False,ax=ax, color='blue', markersize=5)
    rm_t.plot(legend=False,ax=ax, color='black', markersize=5) #Helen wanted non-blue markers
    levee_removal.plot(ax=ax, color='red', edgecolor='red', alpha=0.5)
    
    ctx.add_basemap(ax=ax, source = ctx.providers.Esri.WorldImagery, attribution=False, attribution_size=6,
                    crs = 'epsg:26910', alpha=0.7)
    

    
ax.legend(handles=map_lgd, loc='upper left')

main_map(ax)
plt_arrow(ax)
plt_scale(ax)
mark_inset(axins, ax, loc1=2, loc2=3, fc="none", ec="black")
# fig.savefig(join(fig_dir, 'CA_domain_map.png'), bbox_inches='tight')

# +
fig, ax = plt.subplots(figsize=(6.5, 6.5), dpi=300)

main_map(ax=ax)
fontsize=10
arr_lab(lak_extent, 'Reconnected\nFloodplain', ax, offset = (-100, 150), fontsize=fontsize)
arr_lab(lak_extent, 'Cosumnes\nRiver', ax, offset = (-150, -450), fontsize=fontsize)
arr_lab(m_domain, 'Model\nDomain', ax, offset = (-600, -0), fontsize=fontsize)

ax.legend(handles=lgd_short, loc='upper left')

dir_arrow(ax, 0.8, 0.9, -0.075, -0.075, -0.05, 'Streamflow\nDirection', fontsize=8)

plt_arrow(ax, 0.05, 0.15)
make_multi_scale(ax, 0.1, 0.1, dist = 0.5E3, scales = [4,2,1])

## inset map
axins = inset_axes(ax, width="35%", height="35%", 
                   loc='lower right',
                   # bbox_to_anchor=(-.01, .01, 1, 1),
                   bbox_to_anchor=(-.01, .01, 1, 1),
                  bbox_transform=ax.transAxes, 
                   # loc=2
                  )

ref_map(axins, fontsize=8)
arr_lab(m_domain, 'Cosumnes\nWatershed', axins, offset = (125, 40), arrow=False, fontsize=8)


# -

# ## Interior floodplain plot of wells

# +
fig, ax =plt.subplots(figsize=(6.5,6.5))
ax.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)

plt_bnds = gdf_bnds(rm_t,ax=ax, buf=200)

m_domain.plot(color="none",edgecolor='black',ax=ax)

cr.plot(ax=ax, color='tab:blue')
mr.plot(ax=ax, color='tab:blue')

lak_extent.plot(ax=ax, color='none', hatch='///', edgecolor='white',  alpha=0.7)
rm_t.plot(legend=False,ax=ax, color='black', markersize=5)


rm_t.apply(lambda x: ax.annotate(x.Sensor.replace('MW_',''), xy=x.geometry.coords[0], ha='center', fontsize=6,
                                    xytext = (5,10), textcoords='offset pixels',
                                    bbox=dict(boxstyle="square,pad=0.3", fc="lightgrey", ec="black", lw=2)
                                                        ),axis=1);# ax.annotate?

ctx.add_basemap(ax=ax, source = ctx.providers.Esri.WorldImagery, attribution=False, attribution_size=6,
                crs = 'epsg:26910', alpha=0.7)
# fig.savefig(join(fig_dir, 'floodplain_map.png'), bbox_inches='tight')
# -

# ### Need to identify extent of floodplain on channel

# # Model summary plots
# It would be cool to have functions that are built on using the flopy output rather than using input data even slightly slower because then there is more interoperability

# +
loadpth = 'F:/WRDAPP/GWFlowModel/Cosumnes/Stream_seepage/'
model_nam = 'oneto_denier_upscale4x_2014_2018/'

loadpth = 'C:/WRDAPP/GWFlowModel/Cosumnes/Stream_seepage/'
model_nam = 'oneto_denier_upscale4x_2014_2020_no_reconnection/'

model_ws = loadpth+model_nam
model_ws
# -

m = flopy.modflow.Modflow.load(f = 'MF.nam',  model_ws = model_ws)

dis = m.dis

tmuni = ['s','m','h','D','M','Y'][m.dis.itmuni-1]
tr = ~m.dis.steady.array
delr = m.dis.delr[0]
delc = m.dis.delc[0]


strt_date, end_date, dt_ref = get_dates(m.dis, ref='strt')


# ### Hydrologic input plots (methods)

# +
sfr_dir = join(gwfm_dir, 'SFR_data')
# USGS presents flow in cfs (cubic feet per second)
inflow_in = pd.read_csv(join(sfr_dir,'MB_daily_flow_cfs_2010_2019.csv'), index_col = 'datetime', parse_dates = True)

# covnert flow from cubic feet per second to cubic meters per day
inflow_in['flow_cms'] = inflow_in.flow_cfs * (1/(3.28**3))
# filter out data between the stress period dates
inflow = inflow_in.loc[strt_date:end_date]
# -

## Potential ETo spatial interpolation from CIMIS
fn = glob.glob(join(uzf_dir,'CIMIS', 'Cosumnes_dailyET_precip*.csv'))
daily_data = pd.DataFrame()
for file in fn:
    new_data = pd.read_csv(file, index_col = ['Date'], parse_dates = True)
    daily_data = pd.concat((daily_data, new_data))
# units of mm
data_in = daily_data[daily_data['Stn Name']=='Fair Oaks']
# clean up data so columns are by location, units of Precip are in mm
rain_in = data_in.pivot_table(index = 'Date', columns = 'Stn Name', values = 'Precip (mm)')
rain_m = rain_in/1000
# subset to model and reindex to match streamflow dates
rain_plt = rain_m['Fair Oaks'].loc[strt_date:end_date]
rain_plt = rain_plt.reindex(inflow.index)

flw_lgd =  [
    # Patch(facecolor='none', edgecolor='black', linestyle = '--', alpha=1., label='Watershed Extent'),
    Line2D([0], [0],color='tab:blue', alpha=0.6, label='Rainfall'),
    Line2D([0], [0], color='brown', alpha=0.6,label='Streamflow'),
    # Line2D([0], [0], linestyle='-.', color='black', label='23 cms'),
    # Line2D([0], [0],  linestyle='--', color='black', label='71.6 cms'),
]


# +
fig, ax = plt.subplots(figsize=(6.5,3))

dt = inflow.index.values
ax.plot(dt, inflow['flow_cms'], color='brown', alpha=0.6)
# Create second axes, in order to get the bars from the top you can multiply by -1
ax2 = ax.twinx()
# ax2.bar(dt, -rain_plt.values, 0.9)
ax2.plot(dt, -rain_plt, alpha=0.6)# regular line plot does better than a bar plot
# Now need to fix the axis labels
max_pre = np.max(rain_plt)+0.01
max_pre = 0.12
y2_ticks = np.arange(0, max_pre, 0.02)
y2_ticklabels = [str(i) for i in y2_ticks]
ax2.set_yticks(-1 * y2_ticks)
ax2.set_yticklabels(y2_ticklabels)

ax2.set_ylabel('Rainfall (m)')
ax.set_ylabel('Streamflow ($m^3/s$)')
# ax.ticklabel_format(style='plain', axis='y') 
ax.set_yscale('log')
ax.set_ylim(1,1E3)
ax.set_xlabel('Date')

fig.legend(handles=flw_lgd, loc='outside upper center',ncol=2, bbox_to_anchor=(0.5, 1.0)) # with flow thresholds (0.5, 1.05)
# arr_lab(lak_extent, 'Reconnected\nFloodplain', ax, offset = (-100, 150), fontsize=fontsize)
def flow_threshold_lines(ax):
    ax.axhline(y=23, color='black', linestyle='--')
    ax.axhline(y=71.6, color='black', linestyle='-.')
    xy_lab((pd.to_datetime('2015-8-1'),23),'23 $m^3/s$', offset = (0,5), fontsize=10, bbox=False)
    xy_lab((pd.to_datetime('2015-8-1'),71.6),'71.6 $m^3/s$', offset = (0,5), fontsize=10, bbox=False)
flow_threshold_lines(ax=ax)

fig.savefig(join(fig_dir, 'streamflow_hydrograph_with_rain.png'), bbox_inches='tight')
plt.show()

# -

# ### Precipitation only to show WY type

wyt_sac = pd.read_csv(join(gwfm_dir, 'GHB_data', 'sacramento_WY_types.txt'))
color_dict = {'C':'tab:purple','D':'tab:red','BN':'tab:orange','AN':'yellow','W':'tab:green'}
name_dict = {'C':'Critical','D':'Dry', 'BN':'Below Normal', 'AN':'Above Normal','W':'Wet'}
wyt_sac['color'] = [color_dict[yt] for yt in wyt_sac['Yr-type']]
wyt_sac['name'] = [name_dict[yt] for yt in wyt_sac['Yr-type']]


# wyt for plotting
wy = dt_ref.wy.unique()
wyt = wyt_sac[wyt_sac.WY.isin(wy)].copy()
wyt['plt_date'] = wyt.WY.astype(str)+'-1-1'
wyt = wyt.set_index('WY')

annual_rain = rain_m['Fair Oaks'].resample('AS-Oct').sum()
# look at data from last 20 years
annual_rain = annual_rain[annual_rain.index.year>annual_rain.index.year.max()-20]
# # summarize typical rainfall and devication
print('Avg rain %.3f m' %annual_rain.mean())
annual_rain[annual_rain.index.year.isin([2016,2018])] - annual_rain.mean()

# +
plt_strt = '2014-10-1'
plt_end = '2020-9-30'
fig, ax = plt.subplots(figsize=(6.5,3), dpi=300)
rain_mon = rain_plt.resample('MS').sum()
rain_mon = rain_mon.loc[plt_strt:plt_end]

ax.bar(np.arange(0, len(rain_mon)), rain_mon.values)
ax.set_xticks(np.arange(0, len(rain_mon))[3::12], rain_mon.index[3::12].year);
for n in wy:
    xy_lab((np.where(rain_mon.index==wyt.loc[n, 'plt_date'])[0][0], rain_mon.quantile(0.95))
           ,wyt.loc[n, 'name'].replace(' ','\n'), offset = (0,2), fontsize=10, bbox=True)

ax.set_ylabel('Monthly Rainfall Total (m)');
ax.set_xlabel('Date');
plt.savefig(join(fig_dir, 'monthly_rainfall.png'), bbox_inches='tight')
plt.close()

# +
fig, ax = plt.subplots(figsize=(6.5,3), dpi=300)

inflow_plt = inflow.loc[plt_strt:plt_end]

inflow_plt.plot(y='flow_cms', ax=ax, legend=False)
plt.yscale('log')
plt.ylim(1,1.2E3)
ax.set_ylabel('Daily Streamflow ($m^3/s$)')
ax.set_xlabel('Date')
flow_threshold_lines(ax=ax)

plt.savefig(join(fig_dir, 'daily_streamflow.png'), bbox_inches='tight')
plt.close()
# -

# #### Identify general hydrologic characteristics

print('Days above 23 cms:',(inflow.flow_cms>23).sum())
days_by_wy = (inflow.flow_cms>23).resample('AS-Oct').sum()
days_by_wy[days_by_wy.index.year.isin([2016,2018])]

flow = 71.6
print('Days above '+str(flow)+' cms:',(inflow.flow_cms>flow).sum())
days_by_wy = (inflow.flow_cms>flow).resample('AS-Oct').sum()
days_by_wy[days_by_wy.index.year.isin([2016,2018])]

# ## Time series plots

# In Randy Hanson 2010 paper on FMP applications. ET Groundwater is pink, ET irrig is blue, ET precip is orange. Recharge is green, runoff is brown  
#
# For the water budget figure, it is a little deceving to show the input ET since most of that is not used in the end due to a lack of groundwater.

# ### Water budget (results)

wb_color = pd.read_excel('mf_wb_color_dict.xlsx',sheet_name='owhm_wb_dict', comment='#')
color_dict = wb_color.set_index('flux')['color'].to_dict()
label_dict = wb_color.set_index('flux')['name'].to_dict()

wb, out_cols, in_cols = clean_wb(model_ws, dt_ref)

# +
# import matplotlib.ticker as ticker
# ax[0].set_yscale('log')
# ax[0].yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:,.1E}'.format(x)))
# -

fig,ax= plt.subplots(2,1, sharex=True, figsize=(6.5, 6), dpi=300)
scale = 1E6 # convert to millions
# wb.plot(y='PERCENT_ERROR', ax=ax[0])
wb[out_cols].multiply(1/scale).plot(y=out_cols, ax=ax[-2], legend=True, color=[color_dict[n] for n in out_cols], 
        label=[label_dict[n] for n in out_cols])
wb[in_cols].multiply(1/scale).plot(y=in_cols, ax=ax[-1], legend=True, color= [color_dict[n] for n in in_cols], 
        label=[label_dict[n] for n in in_cols])
plt.xlabel('Date')
fig.supylabel('Flux ($m^3/d$)')
fig.supylabel('Flux (million $m^3/d$)')
fig.tight_layout()
fig.savefig(join(fig_dir, 'water_budget_continuous.png'), bbox_inches='tight')

# ## Map spatial coverage of boundary conditions

gde_row, gde_col = np.where(m.evt.exdp[0].array>2)
gde_grid = grid_p.set_index(['row','column']).loc[list(zip(gde_row+1, gde_col+1))]


# +

def spd_2_arr(sp_data, sp_col, dis):
    """Given the stress_period_data from flopy return the data in an array format
    """
    # convert pumping to array
    arr = np.zeros((dis.nper,dis.nrow,dis.ncol))
    for n in np.arange(0,dis.nper):
        data_n = sp_data[n]
        # only index array if there is data for a stress period
        if data_n is not None:
            arr[n, data_n.i, data_n.j] += data_n[sp_col]
    return(arr)


# -

wel_arr = spd_2_arr(m.wel.stress_period_data, 'flux', m.dis)
wel_rate = wel_arr/(dis.delr[0]*dis.delc[0])
wel_row, wel_col = np.where(wel_rate.sum(axis=0)<0)

# +

wel_ss = pd.DataFrame(m.wel.stress_period_data[0])
ag_grid_p = grid_p.set_index(['row','column']).loc[list(zip(wel_row+1, wel_col+1))].reset_index()
# Helen wants points not polygons
ag_grid_p.geometry = ag_grid_p.geometry.centroid

# save dataframe of stream reach data
# sfrdf = pd.DataFrame(m.sfr.reach_data)
# sfrdf[['row','column']] = sfrdf[['i','j']]+1
# # grid_sfr = grid_match.set_index(['row','column']).loc[list(zip(sfrdf.i+1,sfrdf.j+1))].reset_index()
# grid_sfr = grid_match.merge(sfrdf, on=['row','column']).sort_values('iseg')

grid_sfr = reach_data_gdf(m.sfr, grid_p)
grid_sfr_full = grid_sfr.copy()
# for mapping only keep non-zero hydraulic conductivity reaches (stream-aquifer active)
grid_sfr = grid_sfr[grid_sfr.strhc1!=0]
grid_sfr['plt_seg'] = np.arange(1, len(grid_sfr)+1)



# -


lakrow,lakcol = np.where(m.lak.lakarr[0].array[0]==1)
lak_grid = grid_p.set_index(['row','column']).loc[list(zip(lakrow+1, lakcol+1))].reset_index() # drops interior lines
# lak_grid = grid_match.set_index(['row','column']).loc[list(zip(lakrow+1, lakcol+1))].reset_index() # keeps interior lines of regional grid
lak_grid_union = gpd.GeoDataFrame([0], geometry = [lak_grid.geometry.unary_union], crs=grid_p.crs)


# +
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

legend_elements = [
    Patch(facecolor='brown',alpha=0.8,label='Boundary Groundwater Flow'),
    Patch(facecolor='aqua',alpha=0.8,label='GDE Evapotranspiration'),
    # Patch(facecolor='red', edgecolor='r',alpha=0.6,label='Agricultural Pumping Wells'),
    Line2D([0], [0],color='None', marker='.', markerfacecolor='black', label='Agricultural Pumping Wells', linewidth=4),
    # Patch(facecolor='red', edgecolor='r',alpha=0.6,label='Irrigated Lands'),
    Patch(facecolor='none', edgecolor='white',alpha=0.7, hatch="///", label='Floodplain Cells'),
#     Line2D([0], [0],color='tab:blue',label='Stream Segments', linewidth=4),
    Patch(facecolor='tab:blue', label='Stream Segments'),
    # Line2D([0],[0], color='none', markerfacecolor='black', marker='$10$', label='Stream Segments'),
                    ]


# +
# # make square legend items
# import matplotlib
# matplotlib.rcParams['legend.handlelength'] = 1
# matplotlib.rcParams['legend.handleheight'] = 1
# -

# identify reaches where floodplain connects
fp_in_reach = grid_sfr[grid_sfr.iseg == grid_sfr_full[grid_sfr_full.outseg==-1].iseg.values[0]-2].plt_seg.max()-1
fp_out_reach = grid_sfr[grid_sfr.iseg == grid_sfr_full[grid_sfr_full.iupseg==-1].iseg.values[0]+1].plt_seg.max()-1
fp_in_reach, fp_out_reach

plt_segs = np.append(0,np.arange(9,80,10)) # sfr segments to label
plt_segs = np.sort(np.unique(np.append(plt_segs, (fp_in_reach, fp_out_reach))))
plt_segs
# plt_segs = [0,9,19, fp_in_reach, 39, fp_out_reach, 59,69,79]

# +
fig, ax = plt.subplots(figsize=(6.5,6.5), dpi=300)
# ag_grid_p.plot(ax=ax, color='red', alpha=0.6)
gde_grid.plot(ax=ax, color='aqua', alpha=0.6)
# GDE_grid_p.plot(ax=ax, color='aqua', alpha=0.6)
bnd_cells.plot(ax=ax, color='brown')
m_domain.plot(ax=ax, color='none', edgecolor='black')

ag_grid_p.plot(ax=ax, color='black', markersize=5)

# lak_extent.plot(ax=ax,color='none',edgecolor='white',hatch="///", alpha=0.7)
lak_grid_union.plot(ax=ax,color='none',edgecolor='white',hatch="///", alpha=0.7)
grid_sfr.plot(ax=ax, color='tab:blue')
grid_sfr.iloc[plt_segs].apply(lambda x: ax.annotate(str(x.plt_seg), xy=x.geometry.centroid.coords[0], ha='center', fontsize=8,
                                    xytext = (1,1), textcoords='offset pixels', 
                                    bbox=dict(boxstyle="square,pad=0.1", fc="lightgrey", ec="none", lw=1, alpha=0.5) 
                                               ),axis=1);

# arr_lab(grid_sfr.iloc[[0]], 'Segment \nNumber', ax, offset = (120, 10), arrow=True, fontsize=8)
arr_lab(grid_sfr.iloc[[fp_in_reach]], 'Floodplain\nInflow', ax, offset = (180, -100), arrow=True, fontsize=8)
arr_lab(grid_sfr.iloc[[fp_out_reach]], 'Floodplain\nOutflow', ax, offset = (180, -100), arrow=True, fontsize=8)
arr_lab(grid_sfr.iloc[[9]], 'Reach \nNumber', ax, offset = (100, -180), arrow=True, fontsize=8)

ax.legend(handles=legend_elements, loc='lower right')

ctx.add_basemap(ax=ax, source = ctx.providers.Esri.WorldImagery, attribution=False, attribution_size=6,
                crs = 'epsg:26910', alpha=0.8)

# drop axis labels for cleaner plot
ax.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)

## reference items

# dir_arrow(ax, 0.8, 0.9, -0.075, -0.075, -0.05, 'Streamflow\nDirection', fontsize=8)
plt_arrow(ax, 0.05, 0.15)
make_multi_scale(ax, 0.1, 0.1, dist = 0.5E3, scales = [4,2,1])

# -

# I don't think I ever scaled ET by the area it covers in each cell so it potentially overestimates the ET

# # TPROGs

dem_data = np.loadtxt(gwfm_dir+'\DIS_data\dem_52_9_200m_mean.tsv')


import glob as glob
import tprogs_cleaning as tc
import pyvista as pv


# +

tprogs_id = ''
mf_tprogs_dir = gwfm_dir+'/UPW_data/tprogs_final' + tprogs_id+'/'
tprogs_files = glob.glob(mf_tprogs_dir+'*')
# tprogs_files

# -

t = 0
tprogs_line = np.loadtxt(tprogs_files[t])
# convert any negatives representing input data to same value
tprogs_arr = np.abs(np.reshape(tprogs_line, (320, 100,230)))


# +
# bottom elevations array of tprogs
bot_elev = np.reshape(np.flip(np.arange(-80,80,0.5)), (320, 1,1))
bot_elev = np.repeat(np.repeat(bot_elev, 100, axis=1), 230, axis=2)


tprogs_info = [80, -80, 320]

tprogs_lay = tc.elev_to_tprogs_layers(elev=dem_data, tprogs_info=tprogs_info)

# -

arr_dim = (320, 100, 230)


# Rename pyvista.UniformGrid to pyvista.ImageData, wrap the new class as a shell UniformGrid class and issue a Deprecation warning in its __init__

def mfarr2grid(arr):
    # grid = pv.UniformGrid()
    grid = pv.ImageData()
    # Set the grid dimensions: shape because we want to inject our values on the
    # I have to add 1 to each dimension to have it be built on the cells
    grid.dimensions = [101, 231, 321]
    # real origin, but incorrect because of no rotation
    # simple origin that allows easier data output cleaning
    grid.origin = (0, 0, 0) # bottom left corner of the dataset
    grid.spacing = (200, 200, 0.5)
    arr_in = np.moveaxis(arr, 0, 2).flatten(order='F').astype(int)
    grid.cell_data["facies"] = arr_in

    return(grid)


grid_sfr_p = gpd.read_file(join(gwfm_dir, 'SFR_data', 'final_grid_sfr', 'grid_sfr.shp'))
grid_sfr_p = grid_sfr_p.rename(columns={'row':'row_p','column':'column_p'})


# +
def make_river(grid_sfr, lak_grid, tprogs_lay):
    river_arr = np.zeros(arr_dim)
    r_lay = tprogs_lay[grid_sfr.row_p.astype(int)-1, grid_sfr.column_p.astype(int)-1]
    river_arr[r_lay-4, grid_sfr.row_p.astype(int)-1, grid_sfr.column_p.astype(int)-1] = 1
    
    lak_lay = tprogs_lay[lak_grid.row_p.astype(int)-1, lak_grid.column_p.astype(int)-1]
    river_arr[lak_lay-2, lak_grid.row_p.astype(int)-1, lak_grid.column_p.astype(int)-1] = 1
    
    river_arr = np.flip(river_arr, axis=0)
    river = mfarr2grid(river_arr)
    river = river.threshold(value = [0.9, 1.1], scalars='facies') #, preference='cell'
    return(river)

river = make_river(grid_sfr, lak_grid, tprogs_lay)
# -

# array to multiply others
local_cells = np.zeros(tprogs_arr.shape).astype(bool)
local_cells[:,grid_match.row_p-1, grid_match.column_p-1] = True


# y-rotation: 20 degree is better for regional, 10 for oneto-denier local
def pv_rot(mesh, y_rot = 20):
    mesh.rotate_z(90)
    mesh.rotate_x(10)
    # it seems that the tprogs data is somehow flipped when importing it into pyvista
    # because it requires an extra 180 degree rotation
    mesh.rotate_y(y_rot)
#     mesh.rotate_y(200)


# +

def grid_plt(grid, fig_nam, grid2=None, z_scale = 50, y_rot = 20):
    plotter = pv.Plotter(notebook=False, 
#                          lighting=None,
                         off_screen=True # if true then screenshots work
                        )
#     plotter.show_axes()
#     plotter.show_grid()
#     plotter.show_bounds()
    plotter.background_color='white'
    # show_egdes should be done locally but not regionally
    # but if I add lighting then I might not need edges
    # switched to viridis_r to match hot as high k and cold as low k (gravel is 1 and mud is 4)
    mesh = plotter.add_mesh(grid, scalars="facies", cmap='viridis_r', lighting=False)
    pv_rot(mesh, y_rot = y_rot)
    # 50 x is good for regional, 20x is good for local
    plotter.set_scale(1, 1, z_scale)
    if grid2 is not None:
        mesh = plotter.add_mesh(grid2, color='white')
        pv_rot(mesh, y_rot = y_rot)
    plotter.show(screenshot=fig_nam + '.png')
    


# -

tprogs_figs = join(fig_dir, 'tprogs')

t = 0
# 11 was a realization with good fit
for t in [11]:
# for t in [0,1,2]:
    tprogs_line = np.loadtxt(tprogs_files[t])
    # convert any negatives representing input data to same value
    tprogs_arr = np.abs(np.reshape(tprogs_line, (320, 100,230)))

    tprogs_in = tprogs_arr.copy()
    # crop data above land
    tprogs_in[bot_elev>dem_data] = 0
    # flip to keep orientation for pyvista
    tprogs_in = np.flip(tprogs_in, axis=0)
    # locally plot all facies
    tprogs_local = tprogs_in.copy()
    tprogs_local *= local_cells

    tprogs_local_grid = mfarr2grid(tprogs_local)
    tprogs_local_active = tprogs_local_grid.threshold(value = [0.9, 4.1], scalars='facies') #, preference='cell'

    tprogs_grid = mfarr2grid(tprogs_in)
    tprogs_active = tprogs_grid.threshold(value = [0.9, 4.1], scalars='facies') #, preference='cell'
    grid_plt(tprogs_local_active, join(tprogs_figs,'tprogs_local_facies_r'+str(t).zfill(3)), river, z_scale = 20, y_rot = 10)
    grid_plt(tprogs_active, join(tprogs_figs,'tprogs_facies_r'+str(t).zfill(3)), make_river(grid_sfr_p, lak_grid, tprogs_lay))


# +
grid_plt(tprogs_local_active, join(tprogs_figs,'tprogs_local_facies_r'+str(t).zfill(3)), river, z_scale = 20, y_rot = 10)
grid_plt(tprogs_active, join(tprogs_figs,'tprogs_facies_r'+str(t).zfill(3)), make_river(grid_sfr_p, lak_grid, tprogs_lay))

# on the local scale we can start to distinguish paths but it isn't very clear still
# -


