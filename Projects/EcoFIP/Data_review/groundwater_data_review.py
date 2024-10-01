# +
import sys
import os
from os.path import basename, dirname, join, exists, expanduser
import glob

import pandas as pd
import numpy as np
import geopandas as gpd

import matplotlib.pyplot as plt
import matplotlib as mpl
# -

import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm



# +
usr_dir = expanduser('~')
doc_dir = join(usr_dir,'Documents')
## UCD
gwfm_dir = join(usr_dir,'Box/research_cosumnes/GWFlowModel')

sfr_dir = join(gwfm_dir,'SFR_data')
map_dir = join(gwfm_dir,'Mapping')

## LWA 
db_dir = join(usr_dir,'Dropbox (LWA)', '01_Project-Teams')
proj_dir = join(db_dir, '669.03 - DWR Cosumnes Floodplain Recharge')
fig_dir = join(proj_dir, 'figures')
contour_dir = join(fig_dir, 'head_contours')


# +
git_dir = join(doc_dir,'GitHub','CosumnesRiverRecharge')
def add_path(fxn_dir):
    if fxn_dir not in sys.path:
        sys.path.append(fxn_dir)
        

add_path(join(git_dir, 'python_utilities'))

# -

from gw_contour import export_kde_raster, raster2contours
from map_cln import gdf_bnds, pnt_2_tup, lab_pnt, plt_cln
from report_cln import base_round, fmt

# # Data load
#
#

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

mb_buf = mb_regional.copy()
mb_buf.geometry = mb_buf.geometry.buffer(2000)

cr = rivers_clip[rivers_clip.GNIS_Name=='Cosumnes River']
mr = rivers_clip[rivers_clip.GNIS_Name=='Mokelumne River']

grid_sfr = gpd.read_file(gwfm_dir+'/SFR_data/final_grid_sfr/grid_sfr.shp')


dem_data = np.loadtxt(gwfm_dir+'\DIS_data\dem_52_9_200m_mean.tsv')


plt.imshow(dem_data/0.3048, vmax=250)
# plt.contour(np.flipud(dem_data)/0.3048, levels=np.linspace(np.min(dem_data), np.max(dem_data),10) )
plt.colorbar(shrink=0.5, extend='max', label='GSE (ft amsl)')

# ## DWR Contours
# 1. Manually download historic DWR groundwater contours shapefile
# 2. Crop to project area and save output to project folder
# 3. Delete groundwater contour shapefile but save workflow

# +
# only start in 2012 so they actually don't provide more data than we have in the kriging I did
# contours_all = gpd.read_file(join(usr_dir, 'Downloads', 'i08_GroundwaterDepthSeasonal_Contours', 'i08_GroundwaterDepthSeasonal_Contours.shp' ))

# +
# contours_dom = gpd.overlay(contours_all, mb_buf.to_crs(contours_all.crs))

# +
# contours_dom.MSMT_YEAR.unique()
# -

# ## Periodic GWL data
# It may be necessary to pick out years with the most data to do some contouring or to overlap periods of years with similar conditions (cycles of 3-5 years based on Water Year Type). 

# +
# url = 'https://data.ca.gov/api/3/action/datastore_search?resource_id='
# msmt_id = '231ba173-392a-4f62-91fb-07f4a90746d1'

# +
# import json
# import urllib.request
# msmt_url = url+msmt_id
# fileobj = urllib.request.urlopen(url)
# response_dict = json.loads(fileobj.read())
# # print(response_dict)
# -

ghb_dir = join(gwfm_dir, 'GHB_data')

gwe = pd.read_csv(join(ghb_dir, 'domain_dwr_msmts.csv'))
gwe.msmt_date = pd.to_datetime(gwe.msmt_date, errors='coerce')
# assign year and season for grouping
gwe['year'] = gwe.msmt_date.dt.year
gwe['season'] = 'neither'
gwe.loc[gwe.msmt_date.dt.month.isin([3,4,5]), 'season'] = 'spring'
gwe.loc[gwe.msmt_date.dt.month.isin([9,10,11]), 'season'] = 'fall'

# +
gwe_chk = gwe.groupby(['year','site_code']).count()
# # only keep wells where median year has at least 2 measurements
# keep_sites = gwe_chk[gwe_chk.groupby(['year','site_code']).quantile(0.5).gwe>=2].reset_index().site_code.unique()
# pull wells with more than 20 years of at least 2 msmts
keep_sites = (gwe_chk>=2).groupby('site_code').sum().msmt_date>20
keep_sites = keep_sites[keep_sites].index.values

print(keep_sites.shape, gwe.site_code.unique().shape)
gwe_cln = gwe.loc[gwe.site_code.isin(keep_sites)]
gwe_cln = gwe_cln.dropna(subset='gwe')
# -

# pattern is consistent between fall and spring
gwe_cln.hist(['year'], bins=len(gwe_cln.year.unique())) #.groupby('season')

# +
wyt_sac = pd.read_csv(join(ghb_dir, 'sacramento_WY_types.txt'))
wyt_sac['dry'] = 'tab:blue'
wyt_sac.loc[wyt_sac['Yr-type'].isin(['C','D','BN']), 'dry'] = 'tab:red' 
# wyt_sac
color_dict = {'C':'tab:purple','D':'tab:red','BN':'tab:orange','AN':'yellow','W':'tab:green'}
name_dict = {'C':'Critical','D':'Dry', 'BN':'Below Normal', 'AN':'Above Normal','W':'Wet'}
wyt_sac['color'] = [color_dict[yt] for yt in wyt_sac['Yr-type']]

# calculate the wyt group
wyt_sac['wy_group'] = (wyt_sac['Yr-type'] != wyt_sac['Yr-type'].shift(1)).cumsum()
wyt_sac['dry_group'] = (wyt_sac['dry'] != wyt_sac['dry'].shift(1)).cumsum()
wyt_sac['date'] = pd.to_datetime((wyt_sac.WY-1).astype(str)+'-10-1')
# -




from matplotlib.patches import Patch
wy_lgd = []
for n in color_dict.keys():
    wy_lgd += [Patch(facecolor=color_dict[n], label=name_dict[n], alpha=0.5)]
# wy_lgd
dry_lgd = [Patch(facecolor='tab:red', label='Drier', alpha=0.7),
          Patch(facecolor='tab:blue', label='Wetter', alpha=0.7)]


# +
def plt_wyt(wyt_sac, ax, alpha=0.5):
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    for n in wyt_sac.wy_group.unique():
        wy = wyt_sac.loc[wyt_sac.wy_group==n]
        ax.fill_between([wy.date.min(), wy.date.max()+pd.DateOffset(months=11)], 
                         ylim[1], ylim[0], color=wy.color, alpha=alpha)
    ax.set_xlim(xlim[0],xlim[1])
    ax.legend(handles=wy_lgd, loc='lower right')
    
def plt_dry(wyt_sac, ax, alpha=0.5):
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    for n in wyt_sac.dry_group.unique():
        wy = wyt_sac.loc[wyt_sac.dry_group==n]
        ax.fill_between([wy.date.min(), wy.date.max()+pd.DateOffset(months=11)], 
                         ylim[1], ylim[0], color=wy.dry, alpha=alpha)
    ax.set_xlim(xlim[0],xlim[1])
    # plt.legend(handles=dry_lgd, loc='lower right')


# -

uzf_dir = join(gwfm_dir, 'UZF_data')
# Potential ETo spatial interpolation from CIMIS
fn = glob.glob(join(uzf_dir,'CIMIS','Cosumnes_dailyET_precip*.csv'))
daily_data = pd.DataFrame()
for file in fn:
    new_data = pd.read_csv(file, index_col = ['Date'], parse_dates = True)
    daily_data = pd.concat((daily_data, new_data))
# units of mm
data_in = daily_data[daily_data['Stn Name']=='Fair Oaks']
# clean up data so columns are by location, units of Precip are in mm
rain_in = data_in.pivot_table(index = 'Date', columns = 'Stn Name', values = 'Precip (mm)')
rain_m = rain_in/1000
# create array for every period of rainfall
rain_df = rain_m.resample('D').interpolate('zero')['Fair Oaks']
# only 1997- present

# +
rain_plt = rain_df.copy()

plt_strt = '1999-10-1'
plt_end = '2020-9-30'
fig, ax = plt.subplots(figsize=(6.5,3), dpi=300)
rain_mon = rain_plt.resample('AS-Oct').sum()
rain_mon = rain_mon.loc[plt_strt:plt_end]

ax.bar(np.arange(0, len(rain_mon)), rain_mon.values)
# ax.set_xticks(np.arange(0, len(rain_mon))[3::12], rain_mon.index[3::12].year, rotation=90);
ax.set_xticks(np.arange(0, len(rain_mon))[::], rain_mon.index[::].year+1, rotation=90);
# for n in wy:
#     xy_lab((np.where(rain_mon.index==wyt.loc[n, 'plt_date'])[0][0], rain_mon.quantile(0.95))
#            ,wyt.loc[n, 'name'].replace(' ','\n'), offset = (0,2), fontsize=10, bbox=True)

ax.set_ylabel('Monthly Rainfall Total (m)');
ax.set_xlabel('Date');
# plt.savefig(join(fig_dir, 'monthly_rainfall.png'), bbox_inches='tight')
# plt.close()
# -

strt = pd.to_datetime('1960-10-1')
end = pd.to_datetime('2021-10-1')
# plt_wyt(wyt_sac, (1960,2021))
fig,axes = plt.subplots(2,1,sharex=True)
axes[1].set_xlim(strt,end)
plt_dry(wyt_sac, axes[1])
plt_wyt(wyt_sac, axes[0])


# look for stations near the river with more data to plot in time series
cr_buf = cr.copy()[['GNIS_Name','geometry']]
cr_buf.geometry = cr.buffer(3200)


gwe_gdf = gpd.GeoDataFrame(gwe_cln, geometry= gpd.points_from_xy(gwe_cln.longitude, gwe_cln.latitude), crs='epsg:4326')
gwe_gdf = gwe_gdf.to_crs(grid_sfr.crs)

# wells near river
# 29 within 3200 m, only 9 within 1600 m
# gwe_riv = gpd.overlay(gwe_gdf.to_crs(cr_buf.crs), cr_buf)
# nearest join to track distance
# gwe_riv = gpd.sjoin_nearest(gwe_gdf.to_crs(cr.crs), cr, distance_col = 'riv_dist', max_distance=3200)
gwe_riv = gpd.sjoin_nearest(gwe_gdf.to_crs(grid_sfr.crs), grid_sfr, distance_col = 'riv_dist', max_distance=3200)

gwe_riv_stns = gwe_riv.drop_duplicates('site_code').copy()

gse_quant = gwe_riv_stns.z_min.quantile([0, 0.33,0.66,1]).values
gwe_riv_stns['gwe_group'] = np.nan
for n, gse in enumerate(gse_quant):
    gwe_riv_stns.loc[gwe_riv_stns.z_min>gse, 'gwe_group'] = n+1

gwe_riv_cln = gwe_riv.merge(gwe_riv_stns[['site_code','gwe_group']])

# +
fig,ax = plt.subplots()
grid_sfr.plot(ax=ax)
# cr.plot(ax=ax)
gwe_riv_stns.plot('gwe_group', ax=ax, legend=True, legend_kwds={'label':'Group'})

plt_cln(ax=ax)
# -

import seaborn as sns
sns.relplot(gwe_riv_cln, x='msmt_date', y='gwe',hue='site_code',  
            col='gwe_group', col_wrap=3, legend=False, kind='line', facet_kws={'sharey':True})

long_sites = gwe_plt.groupby('site_code').count()['gwe']
long_sites = long_sites[long_sites>npts].index


gwe_plt = gwe_riv_cln[gwe_riv_cln.msmt_date>'2000-1-1']
npts = 40
long_sites = gwe_plt.groupby('site_code').count()['gwe']
long_sites = long_sites[long_sites>npts].index
gwe_plt = gwe_plt[gwe_plt.site_code.isin(long_sites)]
for n in np.arange(1,4):
    # df_plt = gwe_riv_cln[gwe_riv_cln.gwe_group==n]
    df_plt = gwe_plt[gwe_plt.gwe_group==n]
    fig,ax = plt.subplots()
    for s in df_plt.site_code.unique():
        df_plt[df_plt.site_code==s].plot(x='msmt_date',y='gwe', color='black',
                                         kind='line',ax=ax, legend=False)
    # ax = sns.lineplot(df_plt,x='msmt_date',y='gwe',
    #                   color='black',
    #                   hue='site_code',
    #                   legend=False
    #                  )
    plt_dry(wyt_sac, ax)
    plt.ylabel('Groundwater Elevation (ft amsl)')
    plt.xlabel('Date')
    # plt.savefig(join(fig_dir, 'gw_hydrographs_near_river.png'), bbox_inches='tight')
    plt.show()
    plt.close()

# +
# # grahpic for exit seminar/simplify to a few wells
# fig,ax = plt.subplots(dpi=300)
# ls = ['--','-.',':']
# lc = ['black','gray','lightgray']
# for n in np.arange(1,4):
#     df_plt = gwe_plt[gwe_plt.gwe_group==n]
#     # g_sites = df_plt.site_code.unique()
#     # ordered sites for plotting
#     g_sites = df_plt.groupby('site_code').mean(numeric_only=True).sort_values('gwe').index.values
#     g_sites = g_sites[[int(len(g_sites)/2)]]
#     for s in g_sites:
#         df_plt[df_plt.site_code==s].plot(x='msmt_date',y='gwe', 
#                                          color=lc[n-1], #color='black', linestyle=ls[n-1],
#                                          kind='line',ax=ax, legend=False)

# plt_dry(wyt_sac, ax)
# plt.ylabel('Groundwater Elevation (ft amsl)')
# plt.xlabel('Date')

# -

n=2
# df_plt = gwe_riv_cln[gwe_riv_cln.gwe_group==n]
df_plt = gwe_plt[gwe_plt.gwe_group==n]
# sns.relplot(df_plt, x='msmt_date',y='gwe',
fig,axes=plt.subplots(len(df_plt.site_code.unique()),figsize=(6.5,6.5),dpi=300, sharex=True)
for mw_n, mw in enumerate(df_plt.site_code.unique()):
    ax = axes[mw_n]
    df_plt[df_plt.site_code==mw].plot(x='msmt_date',y='gwe', ax=ax, kind='scatter',
                                      zorder=10, color='black')
    plt_dry(wyt_sac, ax)
    ax.set_ylabel(None)
fig.supylabel('Groundwater Elevation (ft amsl)')
axes[0].set_xlim(df_plt.msmt_date.min(), df_plt.msmt_date.max())
plt.xlabel('Date')
fig.tight_layout(h_pad=-0.1)

# Need to show cycles in groundwater contours, drop by early 1990s and again in 2010s with jumps in early 1980s and 2017.

# +
import pykrige.kriging_tools as kt
from pykrige.ok import OrdinaryKriging
from rasterio.transform import Affine

import rasterio


# -

# define function to interpolate/krige water surface elevations
# part 1/2
def lin_krige(x_in, y_in, z_in, crs,  season,year, folder, res=100, write=True):
    # raster resolution for kriging and output

    gridx = np.arange(np.min(x_in), np.max(x_in), res)
    gridy = np.arange(np.min(y_in), np.max(y_in), res)
    
    # Kriging
    # linear, gaussian, spherical, exponential, hole-effect and power
    OK = OrdinaryKriging(
        x_in,  y_in,  z_in,
        # gaussian overweights low values causing all data to look the same, power looks okay with high lag
        # linear still seems best if I can avoid the singular matrix issue
        variogram_model= 'gaussian', #"linear", 
        # variogram_parameters = {'sill': , 'range': r, 'nugget': 0},
        # variogram_parameters = {'slope': 1E-2, 'nugget': 0},
        verbose=True,
        enable_plotting=True,
    #     exact_values = False,
        enable_statistics = True,
        nlags = 20, # if lags is too low, then higher values seem to dominate?
        pseudo_inv=True
    )

    # z is the kriged grid and ss is the variance grid (sigma ^2)
    z, ss = OK.execute("grid", gridx, gridy)
    # flip data because np sets 0,0 in top left while raster is bottom left
    Z  = np.flip(z.data,axis = 0)
    SS = np.flip(ss.data,axis = 0)

    transform = Affine.translation(np.min(x_in) - res / 2, np.max(y_in) - res / 2) * Affine.scale(res, -res)
    if write == True:
        # Output file creation
        new_dataset1 = rasterio.open(join(ghb_dir, folder, season +str(year)+'_kriged.tif'),
                               mode='w',driver = 'GTiff',
                               height=Z.shape[0], width = Z.shape[1],
                               count = 1, dtype = Z.dtype,
                               crs = crs,  nodata = -9999,  transform = transform)
    
        new_dataset1.write(Z,1)
        new_dataset1.close()
    
        new_dataset2 = rasterio.open(join(ghb_dir, folder, season+str(year)+'_variance.tif'),
                                   mode='w',driver = 'GTiff',
                                   height=Z.shape[0], width = Z.shape[1],
                                   count = 1, dtype = Z.dtype,
                                   crs = crs, nodata = -9999,  transform = transform)
    
        new_dataset2.write(SS,1)
        new_dataset2.close()



season='fall'
year=1965
gwe_in = gwe_gdf[(gwe_gdf.year==year)&(gwe_gdf.season==season)]
gwe_in
ax=gwe_in.plot('gwe')
grid_sfr.plot(ax=ax)


# df = gwe_in.gwe.copy()
def mark_outlier(df):
    """ Series of data"""
    # plot quantiles on the line plot
    quart = df.quantile([.25,.75])
    # median = df.quantile([.5])
    # calculate whiskers
    iqr = quart.loc[0.75]-quart.loc[0.25]
    # 1.5 x the whole interquartile range
    whisker_75 = quart.loc[0.75] + iqr*1.5
    whisker_25 = quart.loc[0.25] - iqr*1.5
    # # where whisker is greater than max or min set as max or min
    whisker_75 =  np.min((df.max(), whisker_75))
    whisker_25 = np.max((df.min(), whisker_25))
    # add column to identify fliers
    df = pd.DataFrame(df).assign(flier= False)
    df.loc[(df.iloc[:,0]<whisker_25)|(df.iloc[:,0]>whisker_75), 'flier']=True
    return(df)
# df = mark_outlier(gwe_in.gwe)
# plt.plot(np.sort(df.loc[df.flier!=True,'gwe'].values), linestyle='',marker='.')


# Removing the outliers caused an odd issue where the groundwater elevations in the foothills were declining toward 0 (average value) further away from the heads with data just below.

wyt_sac[wyt_sac.WY >1960]
wyt_sac.groupby

# +
folder = 'older_interpolated_data'
os.makedirs(join(ghb_dir, folder), exist_ok=True)

season='fall'
for season in ['fall']:
    for year in np.arange(1966, 2022):
        gwe_in = gwe_gdf[(gwe_gdf.year==year)&(gwe_gdf.season==season)]
        #drop outliers
        # gwe_in = gwe_in.loc[mark_outlier(gwe_in.gwe).flier==False]
        # does sorting impact kriging? nope
        # gwe_in = gwe_in.sort_values('gwe').reset_index()
        Z = lin_krige(gwe_in.geometry.x.values, gwe_in.geometry.y.values, gwe_in.gwe.values,
                  gwe_gdf.crs, season, year, folder, write=True, res=100)

# +
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

legend_elements = [Line2D([0], [0], color='b', lw=4, label='Groundwater \nElevation (ft amsl)'),
                   Line2D([0], [0], marker='o', color='black', label='Monitoring Wells',
                          markerfacecolor='black', markersize=8, lw=0),
#                   Patch(facecolor='none', edgecolor='blue',lw=3,
#                          label='Ponds'),
                  # Line2D([0], [0], color='grey', lw=4, label='Groundwater\nFlow Direction' ),
#                   Patch(facecolor='grey', edgecolor='black',lw=3,
#                          label='Groundwater\nFlow Direction'),
                  ]
# -



def plt_contours(gwe_in, gwl_arr, bounds):
    step = 10
    breaks = np.arange(base_round(gwl_arr.min(),base=step)-step, base_round(gwl_arr.max(),base=step)+step,step)

    fig,ax = plt.subplots(figsize=(10.5,8), dpi=300)
    
    # expand map extents and clean up axes
    gdf_bnds(mb_regional, ax=ax, buf=1000)
    # buffer slightly in x direction for legend space
    adjminx, adjmaxx = ax.get_xlim()
    ax.set_xlim(adjminx-1E2, adjmaxx+2E1)
    
    cr.plot(ax=ax, color='tab:blue')
    mb_regional.plot(ax=ax, color='none')
    # plot location of monitoring wells
    # gwe_in.plot(ax=ax, color='red', label='Monitoring Wells')
    gwe_in.plot('gwe', ax=ax, legend=True, vmax=100)
    
    # gwe_in.apply(lambda x: ax.annotate(
    #     # text = str(np.round(x['gwe'], 0)), # head for checking
    #     # text=x['site_code'], # labelling
    #      xy=x.geometry.centroid.coords[0], size=8,
    #     xytext = (-120,-25), textcoords = 'offset pixels', color='black',
    #     bbox=dict(boxstyle="square,pad=0.2", fc="lightgrey", ec="black", lw=2)), axis=1);
    
    ax.legend(handles=legend_elements, loc='lower right')
    
    x_loc, y_loc, arrow_length = 0.05, 0.95, 0.1
    ax.annotate('N', xy=(x_loc, y_loc), xytext=(x_loc, y_loc-arrow_length),
                arrowprops=dict(facecolor='black', width=5, headwidth=15),
                ha='center', va='center', fontsize=20, 
                xycoords=ax.transAxes)
    
    fontprops = fm.FontProperties(size=18)
    scalebar = AnchoredSizeBar(ax.transData,
                               2*5280*0.3048, '2 mi', 'lower left', pad=0.1,sep=5,color='black',
                               frameon=False,size_vertical=200,fontproperties=fontprops)
    ax.add_artist(scalebar)
    plt_cln(ax=ax)
    
    # contour plot options directly from pyplot - needs array and it's extent
    # im = plt.imshow(gwl_arr,alpha=0.5, extent = (bounds.left,bounds.right,bounds.bottom,bounds.top))
    # plt.colorbar(im, label='array')
    # need to flip contours to line up
    CS = plt.contour(np.flipud(gwl_arr), levels= breaks, extent = (bounds.left,bounds.right,bounds.bottom,bounds.top), 
                     colors='blue')
    ax.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=10)
    
    # background imagery Esri.WorldImagery
    # ctx.add_basemap(source=ctx.providers.Esri.WorldImagery, ax=ax, alpha = 0.6,
    #                 crs=lvl_out.crs, attribution=False)
    
    # plt.savefig(join(report_dir,'gw_contour_figure'+str(sdate)+'_2022.png'), bbox_inches='tight')



# +

for season in [ 'fall']:
    for year in np.arange(1966, 2022):
        gwe_in = gwe_gdf[(gwe_gdf.year==year)&(gwe_gdf.season==season)]
        raster_pk = rasterio.open(join(ghb_dir, folder, season +str(year)+'_kriged.tif'))
        gwl_arr = raster_pk.read()[0]
        bounds = raster_pk.bounds
        raster_pk.close()
        plt_contours(gwe_in, gwl_arr, bounds)
        plt.title(str(year))
        plt.savefig(join(contour_dir,'gw_contour_'+season+'_'+str(year)+'.png'), bbox_inches='tight')
        # plt.show()
        plt.close()


# +
from PIL import ImageSequence
from PIL import Image
import glob


for n, season in enumerate(['fall']):
    im_fns = glob.glob(contour_dir+'/*'+season+'*')
    # open first image to start save process then load rest for appending
    im1 = Image.open(im_fns[0]) 
    frames = [Image.open(fn) for fn in im_fns[1:]]
    # write GIF animation
    # takes 30 seconds or so
    fp = open(fig_dir+'/00_contours_'+season+'.gif', "wb")
    # duration is display time of each image in milliseconds, loop = times to repeat
    im1.save(fp, save_all=True, append_images=frames, duration=1000, loop=0)
    fp.close()
# -


