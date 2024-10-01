'''

setback_plotting module
first version: 2024-8-16
'''

import numpy as np
import numpy.ma as ma

import pandas as pd
import geopandas as gpd

import matplotlib.pyplot as plt


def plt_rech(rch_xs_all, flood_type_plot, setback_area=None, log=False, alt_ylab=None):
    rch_sum = pd.DataFrame()
    fig,ax=plt.subplots(3, 1, sharex='col', sharey=False, figsize=(6.5,7), dpi=300)
    
    ft_plt = flood_type_plot.index

    for nf, ft_in in enumerate(ft_plt):
        # second plot for boxplot
        ax_n = ax[nf]    
        rch_xs_sum = rch_xs_all.loc[rch_xs_all.ft==ft_in].drop(columns=['ft'])
        T_in = int(10**flood_type_plot.loc[ft_in,'log_no_d'])

        if setback_area is not None:
            # scale by the number of days to get the daily rate
            rch_xs_sum = rch_xs_sum.multiply(1E6/setback_area, axis=0)/T_in

        rch_xs_sum.transpose().boxplot(ax=ax_n)
        rch_sum = pd.concat((rch_sum, rch_xs_sum.assign(ft=ft_in)))
    for nf, ft_in in enumerate(ft_plt):
        ax_n = ax[nf]
        ax_n.annotate('Type '+str(ft_in),xy=(0.125,0.9),xycoords='axes fraction',
                      bbox={'facecolor': 'lightgray', 'alpha': 0.9, 'pad': 2}
                     )
        if log==True:
            ax_n.set_yscale('log')
    
    # set x labels for boxplots 
    ax_n = ax[-1]
    rot_ticks =  plt.setp(ax_n.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    
    fig.supxlabel('Setback Distance (m)')#,y=-0.04)
    fig.supylabel('Total Recharge (MCM)')#, x=-0.01)
    if setback_area is not None:
        fig.supylabel('Recharge per unit area (m/day)');#, x=-0.01)
    if alt_ylab is not None:
        fig.supylabel(alt_ylab);
    
    fig.tight_layout(h_pad=0.1)
    return rch_sum



def floodplain_annotate(ax, setback_outer, grid_sfr,  title=['(A)', '(B)','(C)']):

    for n in np.arange(0,len(title)):
        ax_n = ax[n]
        ax_n.imshow(setback_outer, cmap='gray') # works when dpi=600
        ax_n.plot(grid_sfr.column-1, grid_sfr.row-1, color='black', linewidth=0.5, linestyle='-.')
        ax_n.annotate(text = title[n], xy=(0.1,0.8), xycoords='axes fraction',
                      bbox={'facecolor': 'lightgray', 'alpha': 0.9, 'pad': 2})

def sfr_setback(grid_sfr, grid_p, s):
    setback = s*200
    nrow,ncol = grid_p[['row','column']].max()
    # grid_sfr = gpd.read_file(gwfm_dir+'/SFR_data/final_grid_sfr/grid_sfr.shp')
    sfr_union = gpd.GeoDataFrame(pd.DataFrame([0]), geometry = [grid_sfr.unary_union], crs=grid_sfr.crs)
    sfr_union.geometry = sfr_union.buffer(setback).exterior
    # grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')
    
    setback_grid = gpd.sjoin(sfr_union, grid_p, how='left')
    setback_outer = np.zeros((nrow,ncol))
    setback_outer[setback_grid.row-1, setback_grid.column-1] = 1
    setback_outer = ma.masked_where(setback_outer==0,setback_outer)
    return(setback_outer)

def arr_plt_label_stream(ax, xy_down=(20,80), xy_up=(180,40), xy_arrow = (140,90)):
    ax.annotate(text = 'Downstream', xy=xy_down, xycoords='data',
                  bbox={'facecolor': 'lightgray', 'alpha': 0.9, 'pad': 2})
    ax.annotate(text = 'Upstream', xy=xy_up, xycoords='data',
                  bbox={'facecolor': 'lightgray', 'alpha': 0.9, 'pad': 2})
    ax.text(xy_arrow[0], xy_arrow[1], s='  N  ',rotation=52.9, size=10,
            bbox=dict(boxstyle="rarrow,pad=0.1",
                      fc="lightgray", ec="black", lw=2))

