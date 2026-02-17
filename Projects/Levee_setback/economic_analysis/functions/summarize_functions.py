# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.6
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import numpy as np
import pandas as pd

# %%
# need to find a better to reference other functions than this
# issue is when calling this function from a different level directory then it fails to find the function
import sys
from os.path import basename, dirname, expanduser, join
usr_dir = expanduser('~')
git_dir = join(usr_dir, 'Documents','GitHub','CosumnesRiverRecharge','Projects','Levee_setback','economic_analysis')
sys.path.append(join(git_dir, 'functions'))
import Basic_soil_budget_monthly as swb


# %%

def format_irr_all(crop, year, crop_in, pred_dict, irr_gw_df_all, irr_sw_df_all):
    """
    Based on a dataframe for gw and sw irrigaton separately (output from model),
    recreate the input that specifies the combined gw and sw irrigation as an array
    """
    var_gen, var_crops, var_yield, season, pred_dict, crop_dict, var_irr = swb.load_var(crop)

    yield_start = swb.ymd2dt(year, season.month_start, season.day_start, season.start_adj)
    yield_end = swb.ymd2dt(year, season.month_end, season.day_end, season.end_adj)
    # get the total extent of the irrigation season (calculation period)
    strt_date = yield_start.min()
    end_date = yield_end.max()
    dates = pd.date_range(strt_date, end_date, freq='D')
    print('Start',strt_date.date(),'end', end_date.date(), 'nper',(end_date-strt_date).days+1)

    gap_irr = var_crops['gap_irr'] # Number of days between irrigations
    n_irr = np.floor(len(dates)/gap_irr).astype(int) + 1 # Calculate number of irrigations
    irr_days = np.arange(0, (n_irr*gap_irr-1), gap_irr).astype(int) # Calculate days on which irrigation takes place
    # specify dates of irrigation to sample
    irr_dates = pd.Series(strt_date+irr_days.astype('timedelta64[D]'))

    print('Preparing irrigation input')
    # get fields for the crop
    crop_df = crop_in[crop_in.name==pred_dict[crop]]
    
    # sample the irrigation rates by crop and irr_date
    irr_gw_crop = irr_gw_df_all.set_index('UniqueID').loc[crop_df.parcel_id].copy()
    irr_gw_crop_dates = irr_gw_crop.reset_index().set_index('date').loc[irr_dates]
    irr_gw_crop_dates = irr_gw_crop_dates.sort_values(['UniqueID','date'])
    irr_gw_crop_dates = irr_gw_crop_dates.reset_index().pivot_table(columns='date', values='rate',index='UniqueID')
    # repeat for Surface water
    irr_sw_crop = irr_sw_df_all.set_index('UniqueID').loc[crop_df.parcel_id].copy()
    irr_sw_crop_dates = irr_sw_crop.reset_index().set_index('date').loc[irr_dates]
    irr_sw_crop_dates = irr_sw_crop_dates.sort_values(['UniqueID','date'])
    # pivot to wide dataframe to match input shape
    irr_sw_crop_dates = irr_sw_crop_dates.reset_index().pivot_table(columns='date', values='rate',index='UniqueID')

    # specify irr_all input
    # irr_gw goes first then irr_sw
    irr_all = np.hstack((irr_sw_crop_dates.values,irr_gw_crop_dates.values))
    return(irr_all)

# %%

def adj_irr_rates(irr_gw_df_all, year, crop_list, irr_est, crop_in, pred_dict, report_change=False):
    """
    Based on data frame of optimized irrigation rates, adjust to account
    for where crop changes require establishment irrigation
    INPUT:
        irr_gw_df_all: dataframe with optimized irrigation rates
        year: current year to reference in crop_in
        crop_list: crops to iterate over
        irr_est: dataframe that defines total irrigation required for the year of establishment
        crop_in: list of parcels with crop and POD information for the given and past years
        pred_dict: translation of crop names
    OUTPUT:
        irr_gw_df_all: dataframe with both optimized irrigation rates and reduced rates where there is establishment
    """
    for crop in crop_list:
        irr_est_c = irr_est[irr_est.crop==crop].copy()
        if len(irr_est_c)>0:
            # find crops that change and need alteration to irr from past year
            for n in irr_est_c.year.unique().astype(int):
                irr_est_n = irr_est_c[irr_est_c.year==n]
                # need to vary this if current or next
                crop_change_n = crop_in.name != crop_in['name_'+(str(year-n))]
                new_crop = crop_in[crop_change_n]
                if n==2:
                    crop_change_n2 = crop_in.name == crop_in['name_'+(str(year-n+1))]
                    new_crop = crop_in[crop_change_n&crop_change_n2]
                if n==3:
                    crop_change_n3 = crop_in.name == crop_in['name_'+(str(year-n+1))]
                    new_crop = crop_in[crop_change_n & crop_change_n2 & crop_change_n3]

                est_irr_id = new_crop[new_crop.name==pred_dict[crop]].parcel_id
                # print(est_irr_id)
                # est_irr_id
                irr_gw_adj = irr_gw_df_all[irr_gw_df_all.UniqueID.isin(est_irr_id)].copy()
                # for now assume pre-harvest and post-harvest are together
                # as the SWB doesn't go past the yield date and the non yield data is from a separate SWB
                # should check back on this
                new_irr_tot = ((irr_est_n.pre_harvest+irr_est_n.post_harvest).values/12)*.3048
                scale_irr_tot = new_irr_tot/irr_gw_adj.groupby('UniqueID')[['rate']].sum()
                # where total irrigation was zero before, make the scalar zero
                scale_irr_tot.loc[scale_irr_tot.rate==np.inf, 'rate'] = 0
                scale_irr_tot = scale_irr_tot.rename(columns={'rate':'irr_adj'})
                irr_gw_adj = irr_gw_adj.merge(scale_irr_tot, on='UniqueID')
                # scale individual irrigation events by the new total
                irr_gw_adj.rate *= irr_gw_adj.irr_adj
                if report_change:
                    # add identifier for years of difference and change
                    irr_gw_adj['changed'] = True
                    irr_gw_adj['year_offset'] = n
                # remove parcels from existing dataframe that are being updated
                irr_gw_df_all = pd.concat((irr_gw_df_all[~irr_gw_df_all.UniqueID.isin(est_irr_id)], irr_gw_adj.drop(columns=['irr_adj'])))
    return(irr_gw_df_all)
