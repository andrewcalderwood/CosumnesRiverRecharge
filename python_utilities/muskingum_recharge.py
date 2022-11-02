"""
muskngum_recharge_functions module. 
Different functions for streamflow routing and stream seepage estimates 
First iteration as a Module July 2022
Author: Andrew Calderwood
"""

import numpy as np
import pandas as pd
import geopandas as gpd 

def Muskingum(I, N, K, X):
    ''' Given an upstream inflow route flow downstream with Muskingum
    I = inflow array 1D
    N = number of reaches
    K = travel time 
    X = wedge storage (range of 0 for none to 0.5 for full wedge, natural stream X is between 0 and 0.3 with a mean of 0.2)
    '''
    T = len(I)
    # N rows and T columns
    dim = [N+1, T]
    dim2 = [N, T]
    Q = np.zeros(dim)
    S = np.zeros(dim2)
    mb = np.zeros(dim2)
    e = np.zeros(dim2)
    
    Q[0,:] = I # set first row as historic inflow
    Q[1:,0] = Q[0,0] # assume initial inflow equals outflow all reaches
    S[:,0] = K*Q[1,0] # assume initial storage equals outflow because no wedge storage
    
    # Muskingum coefficients
    C1 = (delt - 2*K*X)/(2*K*(1-X)+delt)
    C2 = (delt + 2*K*X)/(2*K*(1-X)+delt)
    C3 = (2*K*(1-X)-delt)/(2*K*(1-X)+delt)
#     print('C1:', C1,'C2:', C2, 'C3:', C3)
    
    for j in range(0,T-1):
        for i in np.arange(0, N):
            Q[i+1,j+1] = np.clip(C1*Q[i,j+1] + C2*Q[i,j] + C3*Q[i+1,j], 0, 1E15)
            
    return(Q)

def xs_setback(xs_levee_smooth, setback):
    ''' Function to add levee wall to XS if setback location is not 16 ft elevation above thalweg'''
    mid = np.mean(xs_levee_smooth.index) # location that should be channel bottom based on NHD line
    roll_window = 400 # window used to capture true channel minimum
    # for a given setback imagine there is an impenetrable levee blocking overbank flow
    xs_elevs = xs_levee_smooth[mid-100-setback:mid+100+setback]
    # the channel should fall within the center 400m so that minimum can be used to set new levee height
    thalweg = xs_elevs.loc[mid-roll_window/2:mid+roll_window/2].min()
    # check to see XS elevation at setback distance to determine if it should be raised to needed levee height
    # where XS height is less than 16 ft above channel bottom then raise to 16 ft above
    if xs_elevs.loc[mid-100-setback] - thalweg< 16*0.3048:
        xs_elevs.loc[mid-100-setback] = thalweg + 16*0.3048
    if xs_elevs.loc[mid+100+setback] - thalweg< 16*0.3048:
        xs_elevs.loc[mid+100+setback] = thalweg + 16*0.3048
    return(xs_elevs)

def mannings(d, xs_elevs, n, S):
    ''' Manning Equation from XS data assume constant roughness, slope estimate
    slope = singular value typically based on xs mins
    roughness = singular value, could vary along XS but information is usually not available'''
    # water surface elevation is thalweg plus depth
    wse = xs_elevs.min()+d
    # when wse is above max height of XS then remove excess water (wse = max height)
    if wse > xs_elevs.max():
        wse = xs_elevs.max()
    # calculate depth of water over each XS
    xs_wet = wse - xs_elevs
    # when wse is below XS elevation, remove XS from calculation
    xs_wet[xs_wet < 0] = 0

    # multiply by 10m to get area and sum for area
    A = xs_wet.sum() * 10
    # calculate wetted perimeter
    xy = xs_elevs[wse - xs_elevs > 0]
    # use distance formula to calculate length of XS that is wetted
    Wp = np.sum(np.sqrt(np.diff(xy.values)**2 + np.diff(xy.index.values)**2))
    Q_calc = (np.sqrt(S)/n) * (A**(5/3)) / (Wp**(2/3))
    return(Q_calc)

def min_Q(x, xs_elevs, n, S, Q):
    ''' return error in manning calc to find correct x (depth)'''
    Q_calc = mannings(x, xs_elevs, n, S)
    f0 = abs(Q_calc - Q)
    return(f0)


def mannings_v(x, Q, xs_elevs):
    d = x
    n = 0.048
    S = 0.001  
    wse = xs_elevs.min()+d
    xs_wet = wse - xs_elevs
    xs_wet[xs_wet < 0] = 0
    # multiply by 10m to get area and sum for area
    A = xs_wet.sum() * 10
    # calculate wetted perimeter
    xy = xs_elevs[wse - xs_elevs > 0]
    Wp = np.sum(np.sqrt(np.diff(xy.values)**2 + np.diff(xy.index.values)**2))
    Q_calc = (np.sqrt(S)/n) * (A**(5/3)) / (Wp**(2/3))
    V_calc = Q_calc/A
    return(V_calc)

def gridded_interpolation(vals_in, x_in, y_in, x_out, y_out):
    ''' scipys gridded interpolation
    vals_in are values to be interpolated
    x_in, y_in are locations of values
    x_out, y_out are model grid cell centroids'''
    ncell = nrow*ncol
    # Filter climate data for the strt_date and end_date of the model stress periods
    # Get the xy cell centers for the model grid
    out_xy = np.transpose(np.vstack((x_out, y_out)))

    from scipy.interpolate import griddata
    in_xy = np.transpose(np.vstack((x_in, y_in)))
    # Final reshaped array of interpolated data
    final = np.zeros((nrow, ncol))

    grid = griddata(in_xy, vals_in, xi = out_xy, #in_xy[notna_data]
                    method = 'linear') # maybe try linear?
    # if interpolation fails fill, with other values
    if np.isnan(grid).sum() >0:
        grid_nearest = griddata(out_xy[~np.isnan(grid)], grid[~np.isnan(grid)], xi = out_xy,
                            method = 'nearest')
        grid[np.isnan(grid)] = grid_nearest[np.isnan(grid)]
    final[:,:] = np.reshape(grid, (nrow,ncol))
    return(final)

def calc_depth_arr(wse_grid, wse_all):
    ''' summarize water surface elevation geodataframe to model grid nodes
    wse_grid = geodataframe of join between XS points and model grid
    wse_all is 1D dataframe of the water surface elevation for each XS
    '''
    grid_wse_mean = wse_grid.join(wse_all, on='xs_num')
    grid_wse_mean = grid_wse_mean.dissolve(by = 'node',aggfunc='mean')

    # interpolate flood water surface and elevation from cross-section data
    wse_arr = gridded_interpolation(grid_wse_mean.wse_m, grid_wse_mean.geometry.centroid.x, grid_wse_mean.geometry.centroid.y)
    z_arr = gridded_interpolation(grid_wse_mean.z_m, grid_wse_mean.geometry.centroid.x, grid_wse_mean.geometry.centroid.y)
    # flood depth is difference between flood water surface elevation and cross-section elevation
    d_arr = wse_arr-z_arr
    d_arr[d_arr<0] = np.NaN
    # remove data outside of setback
    d_arr[str_setbacks[s] !=1] = np.NaN
    return(d_arr)    



