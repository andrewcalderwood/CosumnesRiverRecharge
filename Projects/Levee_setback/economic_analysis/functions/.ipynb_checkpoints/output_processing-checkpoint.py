# %%

from os.path import join
import pandas as pd
import numpy as np

import h5py

# %%
# import sys
# sys.path.append('../')


# %%
# import ../functions.swb_functions as swb
# when the module is called it seems to source from a different

# directory than functions
import functions.Basic_soil_budget_monthly as swb

from functions.data_functions import read_crop_arr_h5



# %%
# loadpth = 'C://WRDAPP/GWFlowModel/Cosumnes/Economic'
# m_nam='input_write_2000_2022'
# model_ws = join(loadpth, m_nam)
# year=2010
# # testing
# swb_ws = join(model_ws, 'rep_crop_soilbudget')

# dtw_simple_df = pd.read_csv(join(swb_ws,'field_SWB', 'dtw_ft_WY'+str(year)+'.csv'), index_col='date')
# crop_in = pd.read_csv(join(swb_ws, 'field_SWB', 'crop_parcels_'+str(year)+'.csv'))


# %%
def out_arr_to_long_df(pc_all, out_lin, strt_date, end_date, out_summary):
    """
    Convert the array with data for each field to long format with
    date and UniqueID
    INPUT:
    rch_tot: 2D array with rows of fields and columns of time steps
    out_lin: dataframe that references the represenative output to fields
    strt_date, end_date: are the crop season dates
    """
    # sample represenative data based on dtw_id
    # assumes dtw_id is an integer value that corresponds to the index
    # rch_tot = pc_all[out_lin.dtw_id]
    # # convert pc data to a dataframe with dates
    # pc_df = pd.DataFrame(rch_tot, columns = pd.date_range(strt_date, end_date))

    # need to go from pc_all with simple shape for each representative profile
    # to detailed with infill every 1 ft
    #ach row is a dtw_id and each column is a timestep
    pc_df = pd.DataFrame(pc_all)
    # reindex based on the infilled out_summary
    pc_df = pc_df.reindex(out_summary.dtw_id)
    # linearly interpolate expected rate assuming linearity between every 10 ft dtw
    pc_df = pc_df.interpolate(method='linear')
    # now index using updated index before converting back to array
    pc_df = pc_df.loc[out_lin.dtw_id]
    # convert pc data to a dataframe with dates
    pc_df.columns = pd.date_range(strt_date, end_date)


    # add reference ID columns
    pc_df['UniqueID'] = out_lin.parcel_id
    pc_df['dtw_id'] = out_lin.dtw_id
    # create long format dataframe to append all date to
    pc_df_long = pc_df.melt(id_vars=['UniqueID','dtw_id'], var_name='date',value_name='rate')
    return(pc_df_long)

# %%
# crop = 'Grape'
# var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)
# crop_ref = crop_in[crop_in.name==pred_dict[crop]]


# %%
def downscale_dtw_summary(out_summary):
    """
    given DTW representative values every 10 ft or so,
    provide linear interpolation to give value every 1 foot to provide
    more precise irrigation estimates
    """
    # added if statements to avoid issues probably a cleaner way in the future
    # the downscale goes from 0 to 15 with 0.1 then repeats 15 until end because 16 is missing
    # calculate simple dtw for reindex
    out_summary['round_index'] = out_summary.dtw_mean_ft.round(0)
    if (out_summary.pod_bool==0).any():
        # filter to one type of diversion criteria
        out_summary_gw = out_summary[out_summary.pod_bool==0].set_index('round_index').copy()
        # create new index by going from min to max with smaller step
        out_summary_gw = out_summary_gw.reindex(np.arange(np.min(out_summary_gw.index),np.max(out_summary_gw.index)+1))
        # fill in NA values between 10 ft with linear interpolation
        out_summary_gw = out_summary_gw.interpolate(method='linear').reset_index()
    if (out_summary.pod_bool==1).any():
        out_summary_sw = out_summary[out_summary.pod_bool==1].set_index('round_index').copy()
        out_summary_sw = out_summary_sw.reindex(np.arange(np.min(out_summary_sw.index),np.max(out_summary_sw.index)+1))
        out_summary_sw = out_summary_sw.interpolate(method='linear').reset_index()
    
    if (out_summary.pod_bool==0).any()&(out_summary.pod_bool==1).any():
        out_summary = pd.concat((out_summary_gw, out_summary_sw))
    elif (out_summary.pod_bool==1).any():
        out_summary= out_summary_sw
    elif(out_summary.pod_bool==0).any():
        out_summary = out_summary_gw
    return(out_summary)


    # %%
    # fairly constant across crops
def get_local_data(dtw_simple_df, well_dtw, crop_ref, irr_gw, irr_sw, pc_all):
    """
    Using the representative soil budget output
    transfer the irrigation and percolation rates to individual fields

    OUTPUT:
    out_lin: dataframe that identifies the dtw_id for each UniqueID (parcel)
    out_summary: reference of dtw_id to to dtw_mean_ft
    """
    dtw_df_mean = dtw_simple_df.mean().values
    # if fields constains both with and without pod then double to account alternate run
    if (crop_ref.pod_bool==1).any()&(crop_ref.pod_bool==0).any():
        dtw_df_mean = np.hstack([dtw_df_mean]*2)
        pod_bool = np.repeat([0, 1], int(len(irr_gw)/2))
    elif (crop_ref.pod_bool==1).all():
        pod_bool = 1
    else:
        pod_bool = 0
        
    # # summary output from hdf5 into csv
    out_summary = pd.DataFrame(np.transpose((irr_gw.sum(axis=1), irr_sw.sum(axis=1), pc_all.sum(axis=1), dtw_df_mean)), 
                 columns=['irr_gw_m','irr_sw_m','pc_m', 'dtw_mean_ft'])
    # # actually need to reference the array data directly rather than the sum
    # the nearest value merge should identify the index which should then be used to reference the irrigation/percolation
    # arrays
    out_summary['dtw_id'] = np.arange(0, len(out_summary))
    # specify whether the representative profile is for with or without a POD
    out_summary['pod_bool'] = pod_bool

    # downscale out summary to every 1 ft instead of every 5 ft
    # a simple way would be to reindex based on dtw_mean_ft and infill with a linear range
    out_summary = downscale_dtw_summary(out_summary)
    
    # take the mean dtw for each field
    # dtw_id_mean = dtw_df.mean() # old version continuous
    dtw_id_mean = well_dtw.set_index('UniqueID').dtw_ft # new version ending level at each well
    # sample the dtw for the field simulated for that crop
    dtw_id_mean = pd.DataFrame(dtw_id_mean.loc[crop_ref.parcel_id.values])
    # rename column to match previous coding
    dtw_id_mean = dtw_id_mean.rename(columns={'dtw_ft':'dtw_mean_ft'})
    # want to join dtw info with pod info
    dtw_id_mean = crop_ref.join(dtw_id_mean, on='parcel_id')
    # # add ID to identify which parcel is selected in a simpler id than parcel_id
    dtw_id_mean['id'] = np.arange(0, len(dtw_id_mean))
    
    
    out_summary_dtw = out_summary[['dtw_id','dtw_mean_ft', 'pod_bool']].sort_values('dtw_mean_ft')
    
    # another improvement could be to round dtw_ft to the nearest 0.1 and ID join instead of merge_asof
    
    # identifies the nearest dtw value
    # using the slope in irr/recharge at each point scaled by difference in 
    # need to do separately for pod/no pod
    out_lin0 = pd.merge_asof( 
        dtw_id_mean[dtw_id_mean.pod_bool==0].sort_values('dtw_mean_ft'),
        out_summary_dtw[out_summary_dtw.pod_bool==0].drop(columns='pod_bool'),
                            on='dtw_mean_ft', direction='nearest')
    out_lin1 = pd.merge_asof(
        dtw_id_mean[dtw_id_mean.pod_bool==1].sort_values('dtw_mean_ft'), 
        out_summary_dtw[out_summary_dtw.pod_bool==1].drop(columns='pod_bool'),
                            on='dtw_mean_ft', direction='nearest')
    # rejoin complete table of pod and no pod
    out_lin = pd.concat((out_lin0, out_lin1))
    # sort values for plotting and add back the interpoalted DTW
    out_lin = out_lin.rename(columns={'dtw_mean_ft':'dtw_mean_ft_sim'}).merge( out_summary_dtw)

    return(out_lin, out_summary)

# %%

    # dtw_df_mean = dtw_simple_df.mean().values
    # # if fields constains both with and without pod then double to account alternate run
    # # if (crop_ref.pod_bool==1).any()&(crop_ref.pod_bool==0).any():
    # dtw_df_mean = np.hstack([dtw_df_mean]*2)
    # pod_bool = np.repeat([0, 1], int(len(irr_gw)/2))
    # # elif (crop_ref.pod_bool==1).all():
    # #     pod_bool = 1
    # # else:
    # # pod_bool = 0

    # irr_gw = np.ones(dtw_df_mean.shape)
    # irr_sw = np.zeros(dtw_df_mean.shape)
    # pc_all = np.ones(dtw_df_mean.shape)*0.6
        
    # # # summary output from hdf5 into csv
    # out_summary = pd.DataFrame(np.transpose((irr_gw, irr_sw, pc_all, dtw_df_mean)), 
    #              columns=['irr_gw_m','irr_sw_m','pc_m', 'dtw_mean_ft'])
    # # # actually need to reference the array data directly rather than the sum
    # # the nearest value merge should identify the index which should then be used to reference the irrigation/percolation
    # # arrays
    # out_summary['dtw_id'] = np.arange(0, len(out_summary))
    # # specify whether the representative profile is for with or without a POD
    # out_summary['pod_bool'] = pod_bool
    
    # # take the mean dtw for each field
    # # dtw_id_mean = dtw_df.mean() # old version continuous
    # # dtw_id_mean = well_dtw.set_index('UniqueID').dtw_ft # new version ending level at each well
    # # # sample the dtw for the field simulated for that crop
    # # dtw_id_mean = pd.DataFrame(dtw_id_mean.loc[crop_ref.parcel_id.values])
    # # # rename column to match previous coding
    # # dtw_id_mean = dtw_id_mean.rename(columns={'dtw_ft':'dtw_mean_ft'})
    # # # want to join dtw info with pod info
    # # dtw_id_mean = crop_ref.join(dtw_id_mean, on='parcel_id')
    # # # # add ID to identify which parcel is selected in a simpler id than parcel_id
    # # dtw_id_mean['id'] = np.arange(0, len(dtw_id_mean))
    
    
    # out_summary_dtw = out_summary[['dtw_id','dtw_mean_ft', 'pod_bool']].sort_values('dtw_mean_ft')
    # # a simple way would be to reindex based on dtw_mean_ft and infill with a linear range
    
    # # identifies the nearest dtw value, should correct with a linear interpolation
    # # using the slope in irr/recharge at each point scaled by difference in 
    # # need to do separately for pod/no pod
    # # out_lin0 = pd.merge_asof( 
    # #     dtw_id_mean[dtw_id_mean.pod_bool==0].sort_values('dtw_mean_ft'),
    # #     out_summary_dtw[out_summary_dtw.pod_bool==0].drop(columns='pod_bool'),
    # #                         on='dtw_mean_ft', direction='nearest')
    # # out_lin1 = pd.merge_asof(
    # #     dtw_id_mean[dtw_id_mean.pod_bool==1].sort_values('dtw_mean_ft'), 
    # #     out_summary_dtw[out_summary_dtw.pod_bool==1].drop(columns='pod_bool'),
    # #                         on='dtw_mean_ft', direction='nearest')
    # # # rejoin complete table of pod and no pod
    # # out_lin = pd.concat((out_lin0, out_lin1))
    # # # sort values for plotting and add back the interpoalted DTW
    # # out_lin = out_lin.rename(columns={'dtw_mean_ft':'dtw_mean_ft_sim'}).merge( out_summary_dtw)


# %%
def get_wb_by_parcel(model_ws, year,
                     crop_in, crop_list, dtw_simple_df, well_dtw):
    """ 
    Complete function to take representative soil budget output
    and translate it to each parcel by crop type
    Need to also have a function to include yield and profit
    INPUT:
    
    OUTPUT:
    pc_df_all, irr_sw_df_all, irr_gw_df_all: dataframes with columns for
        UniqueID, date, rate of flux, dtw_id
    """
    # initialize dataframes to save all the output from irrigation and recharge
    pc_df_all = pd.DataFrame()
    irr_sw_df_all = pd.DataFrame()
    irr_gw_df_all = pd.DataFrame()

    for crop in crop_list:
    # for crop in ['Grape']:
        print('Adding', crop,' to the parcel output dataframe')
        var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)

        # need separte hdf5 for each year because total is 300MB, group by crop in array
        fn = join(model_ws, 'field_SWB', "percolation_WY"+str(year)+".hdf5")
        pc_all = read_crop_arr_h5(crop, fn)
        
        # # applied water (GW and SW are separate)
        fn = join(model_ws, 'field_SWB', "GW_applied_water_WY"+str(year)+".hdf5")
        irr_gw = read_crop_arr_h5(crop, fn)
        
        fn = join(model_ws, 'field_SWB', "SW_applied_water_WY"+str(year)+".hdf5")
        irr_sw = read_crop_arr_h5(crop, fn)

        # it would be useful to reference start/end dates for inputting the recharge/pumping
        #  get the dates to update modflow
        yield_start = swb.ymd2dt(year, season.month_start, season.day_start, season.start_adj)
        yield_end = swb.ymd2dt(year, season.month_end, season.day_end, season.end_adj)
        # get the total extent of the irrigation season (calculation period)
        strt_date = yield_start.min()
        end_date = yield_end.max()
        
        crop_ref = crop_in[crop_in.name==pred_dict[crop]]
        
        # return the array that references DTW to the irr and pc for each represetnative
        # soil water budget
        out_lin, out_summary = get_local_data(dtw_simple_df, well_dtw, crop_ref, irr_gw, irr_sw, pc_all)
        
        # create a dataframe with the irr or pc for each UniqueID parcel
        pc_df_long = out_arr_to_long_df(pc_all, out_lin, strt_date, end_date, out_summary)
        irr_gw_df_long = out_arr_to_long_df(irr_gw, out_lin, strt_date, end_date, out_summary)
        irr_sw_df_long = out_arr_to_long_df(irr_sw, out_lin, strt_date, end_date, out_summary)
        
        # # append data for each crop to the total dataframe
        pc_df_all =  pd.concat((pc_df_all, pc_df_long))
        irr_gw_df_all =  pd.concat((irr_gw_df_all, irr_gw_df_long))
        irr_sw_df_all =  pd.concat((irr_sw_df_all, irr_sw_df_long))
    return(pc_df_all, irr_gw_df_all, irr_sw_df_all)

    # %%
    # pc_df_all = pd.DataFrame()


    # # for crop in crop_list:
    # for crop in ['Grape']:
    #     print('Adding', crop,' to the parcel output dataframe')
    #     var_gen, var_crops, var_yield, season, pred_dict, crop_dict = swb.load_var(crop)

    #     # need separte hdf5 for each year because total is 300MB, group by crop in array
    #     fn = join(model_ws, 'rep_crop_soilbudget','field_SWB', "percolation_WY"+str(year)+".hdf5")
    #     pc_all = read_crop_arr_h5(crop, fn)

# %%
# out_summary

# %%
# # out_summary_down = downscale_dtw_summary(out_summary)
# # # out_summary_down[(out_summary_down.dtw_id>14)&(out_summary_down.dtw_id<18)]
# # # out_summary_down
# # # need to go from pc_all with simple shape for each representative profile
# # # to detailed with infill every 1 ft
# # #ach row is a dtw_id and each column is a timestep
# pc_df = pd.DataFrame(pc_all)
# # # reindex based on the infilled out_summary
# pc_df = pc_df.reindex(out_summary.dtw_id)
# # # linearly interpolate expected rate assuming linearity between every 10 ft dtw
# pc_df = pc_df.interpolate(method='linear')
# # now index using updated index before converting back to array
# pc_df.sum().sum()
# # # dtw_id goes uniquely from 0-n_irr or n_irr*2 if pod==0|1
# # out_summary
# # (out_summary.pod_bool==0).any()&(out_summary.pod_bool==1).any()

# # add reference ID columns
# pc_df['UniqueID'] = out_summary.parcel_id
# pc_df['dtw_id'] = out_lin.dtw_id
# # create long format dataframe to append all date to
# pc_df_long = pc_df.melt(id_vars=['UniqueID','dtw_id'], var_name='date',value_name='rate')
