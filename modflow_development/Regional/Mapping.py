# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import os
from os.path import basename, dirname, exists, join
import sys

import geopandas as gpd
import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

import contextily as ctx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

# %%
## Set up directory referencing
# Package data
usr_dir = os.getcwd()
while basename(usr_dir)!='Users':
    temp = basename(usr_dir)
    usr_dir = dirname(usr_dir)
usr_dir = join(usr_dir, temp)
gwfm_dir = join(usr_dir,'Box/research_cosumnes/GWFlowModel')
doc_dir = join(usr_dir, 'Documents')

map_dir = join(gwfm_dir,'Mapping')
sfr_dir = join(gwfm_dir,'SFR_data')


plt_dir = join(map_dir,'figures')

# %%
from importlib import reload


# %%
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
# flopy github path - edited
add_path(doc_dir+'/GitHub/flopy')
import flopy 

# other functions
py_dir = join(doc_dir,'GitHub/CosumnesRiverRecharge/python_utilities')
add_path(py_dir)

from mf_utility import get_layer_from_elev
from map_cln import gdf_bnds, plt_cln, make_multi_scale

# %%
# reload(map_cln)

# %%
grid_sfr = gpd.read_file(join(sfr_dir, 'final_grid_sfr/grid_sfr.shp'))


# %%
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')


# %%
m_domain = gpd.read_file(gwfm_dir+"/DIS_data/NewModelDomain/GWModelDomain_52_9deg_UTM10N_WGS84.shp")

rivers = gpd.read_file(join(sfr_dir,"NHD_MajorRivers/NHD_MajorRivers.shp")).to_crs(m_domain.crs)
rivers_clip = gpd.clip(rivers, m_domain)

soam = gpd.read_file(map_dir+"/so_am_subbasin/so_am_subbasin.shp").to_crs(m_domain.crs)
cos = gpd.read_file(map_dir+"/cos_subbasin/cos_subbasin.shp").to_crs(m_domain.crs)
ca = gpd.read_file(map_dir+"/ca_state_boundary/CA_State_TIGER2016.shp").to_crs(m_domain.crs)


# %%
def regional_scale_arrow(ax):
    x, y, arrow_length = 0.7, 0.15, 0.1
    ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
                arrowprops=dict(facecolor='black', width=5, headwidth=15),
                ha='center', va='center', fontsize=20,
                xycoords=ax.transAxes)

    ax.ticklabel_format(style='plain')


    fontprops = fm.FontProperties(size=18)
    scalebar = AnchoredSizeBar(ax.transData,
                               10000, '10 km', 'lower right', 
                               pad=0.1,sep=2, color='black',
                               frameon=False,size_vertical=5E2, fontproperties=fontprops)
    ax.add_artist(scalebar)


# %%
def regional_arrow(ax, xoff, yoff):
    x, y, arrow_length = xoff, yoff, 0.1
    ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
                arrowprops=dict(facecolor='black', width=5, headwidth=15),
                ha='center', va='center', fontsize=20,
                xycoords=ax.transAxes)
    

fig,ax = plt.subplots(figsize=(6.5,6.5), dpi=300)

m_domain.plot(ax=ax,color="none",edgecolor='black')
gdf_bnds(m_domain,ax)

rivers.loc[rivers.GNIS_Name.isin(['Mokelumne River','South Mokelumne River'])].plot(ax=ax)
cr = rivers.loc[rivers.GNIS_Name=='Cosumnes River']
cr.plot( ax=ax,label='Cosumnes River')

ax.annotate(text='Cosumnes\n River', xy=list(cr.geometry.iloc[10].centroid.coords)[0], 
            xytext = (6,6), textcoords = 'offset pixels',
            bbox=dict(boxstyle="square,pad=0.3", fc="lightgrey", ec="black", lw=2))


ctx.add_basemap(ax, source = ctx.providers.Esri.WorldImagery, crs='epsg:26910', alpha = 0.8, attribution=False)

# regional_scale_arrow(ax=ax)
regional_arrow(ax, 0.65, 0.15)
make_multi_scale(ax, 0.75,0.1, dist=2E3)
plt_cln(ax=ax)

# first one is CA map
axins = inset_axes(ax, width="80%", height="80%", bbox_to_anchor=(0, .5, .5, .5),
                  bbox_transform=ax.transAxes, loc=2)
# second one is smaller inset
axins2 = inset_axes(ax, width="40%", height="40%", bbox_to_anchor=(0.18, .5, .5, .5),
                  bbox_transform=ax.transAxes, loc=2)
axins.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)
axins2.tick_params(labelleft=False, labelbottom=False, left = False, bottom = False)

ca.plot(ax = axins,alpha = 0.2)
soam.plot(ax = axins, edgecolor = 'black')
cos.plot(ax = axins, edgecolor = 'black')

cos.plot(ax = axins2, edgecolor = 'black')
soam.plot(ax = axins2, edgecolor = 'black')
m_domain.plot(ax = axins2, edgecolor = 'black', color = 'none')


# plt.savefig(join(fig_dir, 'regional_domain_map.png'),  bbox_inches='tight')



# %% [markdown]
# # Plot Major Land Uses (GW Pumping)

# %%
uzf_dir = join(gwfm_dir, 'UZF_data')

lu_ag = gpd.read_file(join(uzf_dir, 'county_landuse', 'domain_ag_lu_2018.shp'))
lu_ag = gpd.overlay(lu_ag, m_domain)
lu_ag['area_m2'] = lu_ag.geometry.area
wel_dir = join(gwfm_dir, 'WEL_data')

# load data of locations of domestic wells
dom_loc = pd.read_csv(join(wel_dir, 'ag_res_parcel_domestic_wells.csv'), index_col=0)
# make row,column 0 based
dom_loc.row = (dom_loc.row-1).astype(int)
dom_loc.column = (dom_loc.column -1).astype(int)
# aggregate to the cell level, summing area will keep water usage scaling correct
dom_loc = dom_loc.groupby(['node','row','column', 'CITY']).sum(numeric_only=True).reset_index()
# join to grid for plotting
dom_grid = grid_p.merge(dom_loc,on=['row','column'])

# %%
# not super helpful since dometsic well parcel mapping will cover most
# lu_urban = gpd.read_file(join(uzf_dir, 'county_landuse', 'domain_urban_lu_2018.shp'))


# %%
lu_crop_sum = lu_ag.groupby('name').sum(numeric_only=True).reset_index()

min_area = 1000*43560*0.3048**2 # minimum area to include in plot is 1000 acres
min_area = 200*200*100*230*0.01 # minimum area is 1% of the domain
print('Min area %.2e m^2' %min_area,'vs number of cells %.i' %(min_area/(200**2)), 
      'vs percent of domain %.2f' %(min_area*100/(200**2)/(100*230) ))
main_crops = lu_crop_sum[lu_crop_sum.area_m2 > min_area].name.unique()

# %%
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib import cm


# %%
fig,ax = plt.subplots(figsize=(6.5,6.5), dpi=300)

#normalize item number values to colormap
norm = mpl.colors.Normalize(vmin=0, vmax=len(main_crops))

m_domain.plot(ax=ax,color="none",edgecolor='black')
gdf_bnds(m_domain,ax)


custom_lgd = []
for nc, c in enumerate(main_crops):
    lu_ag[lu_ag.name.isin([c])].plot(ax=ax, color=cm.viridis(norm(nc)), label=c)
    custom_lgd += [Patch(facecolor=cm.viridis(norm(nc)), edgecolor='None', label=c)]

# rural residences
nc +=1
dom_grid.plot(ax=ax, label='Rural residences', color=cm.viridis(norm(nc)))
custom_lgd += [Patch(facecolor=cm.viridis(norm(nc)), edgecolor='None', label='Rural residences')]
# stream grid cells
grid_sfr.plot(ax=ax)
custom_lgd += [Line2D([0], [0], color='tab:blue', lw=4, label='Stream Cells')]

# regional_scale_arrow(ax=ax)
regional_arrow(ax, 0.65, 0.15)
make_multi_scale(ax, 0.75,0.1, dist=2E3)
plt_cln(ax=ax)

ax.legend(handles=custom_lgd, facecolor='lightgray', framealpha=0.6, loc='upper left')

plt.savefig(join(fig_dir, 'regional_major_land_use_map.png'),  bbox_inches='tight')


# %% [markdown]
# ## Map spatial coverage of boundary conditions

# %%
loadpth= 'C:/WRDAPP/GWFlowModel/Cosumnes/Regional'
model_nam = 'Historical_simple_geology_reconnection'
model_ws = join(loadpth, model_nam)

load_only=['DIS','OC','WEL','EVT', 'SFR']
m = flopy.modflow.Modflow.load(model_ws+'/MF.nam', load_only=load_only)

# %%
# get exterior polyline of model grid
grid_bnd = gpd.GeoDataFrame(pd.DataFrame([0]), geometry = [grid_p.unary_union.exterior], crs=grid_p.crs)
# find cells that construct the model boundary
bnd_cells_df = gpd.sjoin(grid_p, grid_bnd)
bnd_cells = bnd_cells_df[['row','column']] - 1
bnd_cells['grid_id'] = np.arange(0,len(bnd_cells))
bnd_rows, bnd_cols = bnd_cells.row.values, bnd_cells.column.values

# %%
lak_shp = join(gwfm_dir,'LAK_data/floodplain_delineation')
lak_grid = gpd.read_file(join(lak_shp, 'lak_grid_cln.shp'))


# %%
dis= m.dis

# %%
# map locations with EVT
evtr_ss = m.evt.evtr.array[0,0]
gde_grid = grid_p.copy()
gde_grid['evtr'] = evtr_ss[gde_grid.row-1, gde_grid.column-1]
gde_grid = gde_grid[gde_grid.evtr!=0]


# %%

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

wel_arr = spd_2_arr(m.wel.stress_period_data, 'flux', m.dis)
wel_rate = wel_arr/(dis.delr[0]*dis.delc[0])
wel_row, wel_col = np.where(wel_rate.sum(axis=0)<0)

# %%

wel_ss = pd.DataFrame(m.wel.stress_period_data[0])
ag_grid_p = grid_p.set_index(['row','column']).loc[list(zip(wel_row+1, wel_col+1))].reset_index()

# save dataframe of stream reach data
sfrdf = pd.DataFrame(m.sfr.reach_data)
grid_sfr = grid_p.set_index(['row','column']).loc[list(zip(sfrdf.i+1,sfrdf.j+1))].reset_index()



# %%
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

legend_elements = [
    Patch(facecolor='brown',alpha=0.8,label='Boundary Groundwater Flow'),
    Patch(facecolor='aqua',alpha=0.8,label='GW ET Possible'),
    Patch(facecolor='red', edgecolor='r',alpha=0.6,label='Pumping Wells'),
    # Patch(facecolor='red', edgecolor='r',alpha=0.6,label='Irrigated Lands'),
    Patch(facecolor='tab:green',alpha=1,label='Reconnected Floodplain'),
#     Line2D([0], [0],color='tab:blue',label='Stream Segments', linewidth=4),
    Patch(facecolor='tab:blue', label='Stream Segments'),
                    ]


# %%
# m_domain.crs

# %%
fig, ax = plt.subplots(figsize=(6.5,6.5), dpi=300)
ag_grid_p.plot(ax=ax, color='red', alpha=0.6)
gde_grid.plot(ax=ax, color='aqua', alpha=0.6)
bnd_cells_df.plot(ax=ax, color='brown')
m_domain.plot(ax=ax, color='none', edgecolor='black', linewidth=0.5)

grid_sfr.plot(ax=ax, color='tab:blue')
lak_grid.plot(ax=ax,color='tab:green')

ax.legend(handles=legend_elements, loc='upper left')

ctx.add_basemap(ax=ax, source = ctx.providers.Esri.WorldImagery, attribution=False, attribution_size=6,
                crs = 'epsg:26910', alpha=0.8)

# drop axis labels for cleaner plot
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

regional_arrow(ax, 0.65, 0.15)
make_multi_scale(ax, 0.75,0.1, dist=2E3)
# regional_scale_arrow(ax)

plt.savefig(join(fig_dir, 'model_boundary_conditions.png'), bbox_inches='tight')

# %% [markdown]
# # TPROGs

# %%
add_path(join(doc_dir, 'GitHub','CosumnesRiverRecharge', 'tprogs_utilities'))

# %%
dem_data = np.loadtxt(gwfm_dir+'\DIS_data\dem_52_9_200m_mean.tsv')


import glob as glob
import tprogs_cleaning as tc
import pyvista as pv


# %%

tprogs_id = ''
mf_tprogs_dir = gwfm_dir+'/UPW_data/tprogs_final' + tprogs_id+'/'
tprogs_files = glob.glob(mf_tprogs_dir+'*')
# tprogs_files

t = 0
tprogs_line = np.loadtxt(tprogs_files[t])
# convert any negatives representing input data to same value
tprogs_arr = np.abs(np.reshape(tprogs_line, (320, 100,230)))


# %%
# bottom elevations array of tprogs
bot_elev = np.reshape(np.flip(np.arange(-80,80,0.5)), (320, 1,1))
bot_elev = np.repeat(np.repeat(bot_elev, 100, axis=1), 230, axis=2)


tprogs_info = [80, -80, 320]

tprogs_lay = tc.elev_to_tprogs_layers(elev=dem_data, tprogs_info=tprogs_info)

arr_dim = (320, 100, 230)

def mfarr2grid(arr):
    grid = pv.UniformGrid()
    # Set the grid dimensions: shape because we want to inject our values on the
    # I have to add 1 to each dimension to have it be built on the cells
    grid.dimensions = [101, 231, 321]
    # real origin, but incorrect because of no rotation
    # simple origin that allows easier data output cleaning
    grid.origin = (0, 0, 0) # bottom left corner of the dataset
    grid.spacing = (200,200,0.5)
    arr_in = np.moveaxis(arr,0,2).flatten(order='F').astype(int)
    grid.cell_data["facies"] = arr_in

    return(grid)



# %%
river_arr = np.zeros(arr_dim)
r_lay = tprogs_lay[grid_sfr.row.astype(int)-1, grid_sfr.column.astype(int)-1]
river_arr[r_lay-2, grid_sfr.row.astype(int)-1, grid_sfr.column.astype(int)-1] = 1
lak_lay = tprogs_lay[lak_grid.row.astype(int)-1, lak_grid.column.astype(int)-1]
river_arr[lak_lay-2, lak_grid.row.astype(int)-1, lak_grid.column.astype(int)-1] = 1

river_arr = np.flip(np.flip(river_arr, axis=0), axis=1)
river = mfarr2grid(river_arr)
river = river.threshold(value = [0.9, 1.1], scalars='facies') #, preference='cell'

# array to multiply others
# local_cells = np.zeros(tprogs_arr.shape).astype(bool)
# local_cells[:,grid_match.row-1, grid_match.column-1] = True


def pv_rot(mesh):
    mesh.rotate_z(90)
    mesh.rotate_x(10)
    # it seems that the tprogs data is somehow flipped when importing it into pyvista
    # because it requires an extra 180 degree rotation
    mesh.rotate_y(10)
#     mesh.rotate_y(200)


# %%
tprogs_figs = join(fig_dir, 'tprogs')
os.makedirs(tprogs_figs, exist_ok=True)


# %%

def grid_plt(grid, fig_nam, grid2=None):
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
    mesh = plotter.add_mesh(grid, scalars="facies", cmap='viridis', lighting=True)
    pv_rot(mesh)
    # 50 x is good for regional, 20x is good for local
    plotter.set_scale(1, 1, 20)
    if grid2 is not None:
        mesh = plotter.add_mesh(grid2, color='black')
        pv_rot(mesh)
    plotter.show(screenshot=fig_nam + '.png')
    

t = 0
# 11 was a realization with good fit
for t in [11]:#[0,1,2]:
    tprogs_line = np.loadtxt(tprogs_files[t])
    # convert any negatives representing input data to same value
    tprogs_arr = np.abs(np.reshape(tprogs_line, (320, 100,230)))

    tprogs_in = tprogs_arr.copy()
    # crop data above land
    tprogs_in[bot_elev>dem_data] = 0
    # flip to keep orientation for pyvista
    tprogs_in = np.flip(tprogs_in, axis=0)

    tprogs_grid = mfarr2grid(tprogs_in)
    tprogs_active = tprogs_grid.threshold(value = [0.9, 4.1], scalars='facies') #, preference='cell'
    grid_plt(tprogs_active, join(tprogs_figs,'tprogs_facies_r'+str(t).zfill(3)), river)


# %%
grid_plt(tprogs_active, join(tprogs_figs,'tprogs_facies_r'+str(t).zfill(3)), river)


