"""
tprogs_cleaning module. 
Different functions for preparing data for use in MODFLOW and other groundwater modeling software.
First iteration as a Module Oct 2021
Author: Andrew Calderwood
"""

import sys
import numpy as np

def tprogs_cut_elev(tprogs_line, dem_data, **kwargs):
    """
    Parameters
    ----------
    tprogs_line : output from TPROGs of line data formatted to be converted by setting z then x then y
    dem_data : 2D array of elevation data of ground surface above which TPROGs should not be real
    rows : number of rows in the TPROGs model
    cols : number of columns in the TPROGs model
    """
    rows = kwargs.get('rows', np.where(np.ones(dem_data.shape)==1)[0])
    cols = kwargs.get('cols', np.where(np.ones(dem_data.shape)==1)[1])
    tprogs_arr = np.reshape(tprogs_line, (320, 100,230))
    tprogs_c = np.reshape(tprogs_arr[:, rows,cols],
                             (tprogs_arr.shape[0],dem_data.shape[0],dem_data.shape[1]))
    tprogs_elev = np.copy(tprogs_c)
    # the bottom layer of the tprogs model is at -50 m amsl and the top layer is 50 m amsl
    t = 0
    for k in np.arange(-80,80,0.5):
        tprogs_elev[t,dem_data<k]= np.NaN
        t+=1
    masked_tprogs = ma.masked_invalid(tprogs_elev)
    return(masked_tprogs)


def int_to_param(tprogs, params):
    """
    Parameters
    ----------
    tprogs: 3D masked array of TPROGs realization
    params: Reference table connecting TPROGs facie to a hydraulic value
    """
    tprogs[tprogs<0] *= -1
    tprogs = tprogs.astype(float)
    # flip tprogs model along z axis to match modflow definition of 0 as top (TPROGS says 0 is bottom)
    tprogs = np.flip(tprogs,axis=0)
    tprogs_K = np.copy(tprogs)
    tprogs_Sy = np.copy(tprogs)
    tprogs_Ss = np.copy(tprogs)
    # hydraulic parameters from fleckenstein 2006
    # I-IV gravel, sand, muddy sand, mud
    # K in m/s, Sy, Ss
    for n in np.arange(1,5):
        tprogs_K[tprogs==n]= params.loc[n,'K_m_d']
    for n in np.arange(1,5):
        tprogs_Sy[tprogs==n]= params.loc[n,'Sy']
    for n in np.arange(1,5):
        tprogs_Ss[tprogs==n]= params.loc[n,'Ss']
            
    return(tprogs_K,tprogs_Sy,tprogs_Ss)


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
    elev_indices = top/lay_thick - elev_round*(1/lay_thick) 
    return(elev_indices.astype(int))


def get_tprogs_for_elev(tprogs_arr, top_elev, bot_elev, **kwargs):
    """
    Function to grab the TPROGs layers by elevation filters and returns
    a 3D array with uneven numbers of filled values in the vertical direction.

    Parameters
    ----------
    tprogs_arr : 3D masked array of TPROGs_realziation
    top_elev : 2D array of elevation setting top reference point
    bot_elev: 2D array of elevation setting bottom reference poing
    """
    rows = kwargs.get('rows', np.where(np.ones(top_elev.shape)==1)[0])
    cols = kwargs.get('cols', np.where(np.ones(top_elev.shape)==1)[1])
    top_indices = elev_to_tprogs_layers(top_elev)
    bot_indices = elev_to_tprogs_layers(bot_elev)
    # find tprogs layer for desired rows and columns
    top_indices = top_indices[rows, cols].astype(int)
    bot_indices = bot_indices[rows, cols].astype(int)
    # the first row of the array will be the top layer and will progress downward until the max bottom is reached
    # with NaNs for rows,cols where there are less layers indexed than the max
    tprogs_subset = np.full(shape = (np.max(bot_indices - top_indices).astype(int), len(rows)),
                       fill_value = np.nan, dtype = float)
    max_layers = np.max(bot_indices - top_indices)
    for k in np.arange(0,max_layers):
        layexist = (bot_indices-top_indices) > k # pick where data should be referenced
        tprogs_subset[k, layexist] = K[top_indices[layexist]+k, rows[layexist], cols[layexist]]
    # return grabbed data in array format if entire domain was used
    if len(rows) == top_elev.shape[0]*top_elev.shape[1]:
        tprogs_subset = np.reshape(tprogs_subset, (max_layers, top_elev.shape[0], top_elev.shape[1]))
    return(tprogs_subset)