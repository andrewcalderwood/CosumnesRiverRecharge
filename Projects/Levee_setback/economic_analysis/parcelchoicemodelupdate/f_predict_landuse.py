# -*- coding: utf-8 -*-
# %%
"""
Created on Wed Nov  8 16:52:58 2023

@author: cameron.speir
Updated to function 1/23/2024
"""


import pandas as pd
import numpy as np
import os

# %%
# Function to get a specific coefficient
def pull_coef(df, cropname, varname):
    coef = df.loc[df['Crop_Eq'] == cropname, varname].values[0]
    return coef

# %%
def predict_crops(data, rev_prior_yr_df, logit_coefs, return_prob=False):
    ''' Predict crop choice based on input variables
    data: dataframe with variables as columns and rows for each parcel
        'parcel_id', 'year', 'crop_cat', 'acres', 'elev_mean',
       'slope_perc_mean', 'dtwfa', 'dtwsp', 'pod'
    rev_prior_yr_df: dataframe of revenue expected for each crop
    logit_coefs: dataframe with the logit coefficients for each variable in data for each crop
    '''
    # Get a vector of crop equations
    crops = logit_coefs['Crop_Eq'].tolist()

    # Add water year critical or dry indicator
    data['dtwfa2'] = data['dtwfa'] * data['dtwfa']
    data['i_pod_wy'] = data['pod'] * data['wy_dc']
    
    # Get the coefficient on prior year revenue - same across crops
    coef_rev_prior_yr = logit_coefs.loc[0, 'rev_prior_yr']
    
    # Revenue for the previous year by crop
    rpy_vector = rev_prior_yr_df.set_index('Crop_Eq')['rev_prior_yr'].to_dict()
    
    # Calculate predicted probabilities
    # pasture is the base category that fills the probability toward 1
    data['N.Mixed_pasture'] = 1
    for crop in crops:
        data[f'N.{crop}'] = np.exp(
            coef_rev_prior_yr * rpy_vector[crop] +
            pull_coef(logit_coefs, crop, 'intercept') +
            pull_coef(logit_coefs, crop, 'acres') * data['acres'] +
            pull_coef(logit_coefs, crop, 'elev_mean') * data['elev_mean'] +
            pull_coef(logit_coefs, crop, 'slope_perc_mean') * data['slope_perc_mean'] +
            pull_coef(logit_coefs, crop, 'pod') * data['pod'] +
            pull_coef(logit_coefs, crop, 'wy_dc') * data['wy_dc'] +
            pull_coef(logit_coefs, crop, 'i_pod_wy') * (data['pod'] * data['wy_dc']) +
            pull_coef(logit_coefs, crop, 'dtwfa') * data['dtwfa'] +
            pull_coef(logit_coefs, crop, 'dtwfa2') * (data['dtwfa']**2)
        )
    
    # need to make sure crop selection follows the same column order used for cumulative probability
    crops_N = data.loc[:,data.columns.str.contains('N.')].columns.str.replace('N.','').values
    # crops_N = logit_coefs['Crop_Eq'].tolist()

    # Calculate Denominator
    data['Denominator'] = 1 + data[[f'N.{crop}' for crop in crops_N]].sum(axis=1)
        
    # Calculate predicted probabilities
    for crop in crops_N:
        data[f'PP.{crop}'] = data[f'N.{crop}'] / data['Denominator']

    if return_prob:
        return data
        
    # Find the maximum probability
    data['maxprob'] = data[[f'PP.{crop}' for crop in crops_N]].max(axis=1)
    
    # Assign Crop Choice based on maximum probability
    conditions = [data[f'PP.{crop}'] == data['maxprob'] for crop in crops_N]
    data['Crop_Choice'] = np.select(conditions, crops_N, default='Unclassified_fallow')
    # a suggestion by the peer review group was that a farmer may not always
    # select the highest probability, thus a random number should be used to select
    # from the distribution to determine the crop used
    # wait for Yusuke to run test with validation before updating
    
    # MODFLOW/irrigation choice info: parcel_id and Crop_Choice
    # not going to extract probabilities because that should be reviewed/processed elsewhere
    return data[['parcel_id', 'Crop_Choice']]

# %%

def predict_crop_probs(data, rev_prior_yr_df, logit_coefs):
    ''' Predict crop choice based on input variables
    data: dataframe with variables as columns and rows for each parcel
        'parcel_id', 'year', 'crop_cat', 'acres', 'elev_mean',
       'slope_perc_mean', 'dtwfa', 'dtwsp', 'pod'
    rev_prior_yr_df: dataframe of revenue expected for each crop
    logit_coefs: dataframe with the logit coefficients for each variable in data for each crop
    '''
    # Get a vector of crop equations
    crops = logit_coefs['Crop_Eq'].tolist()

    # Add water year critical or dry indicator
    data['dtwfa2'] = data['dtwfa'] * data['dtwfa']
    data['i_pod_wy'] = data['pod'] * data['wy_dc']
    
    # Get the coefficient on prior year revenue - same across crops
    coef_rev_prior_yr = logit_coefs.loc[0, 'rev_prior_yr']
    
    # Revenue for the previous year by crop
    rpy_vector = rev_prior_yr_df.set_index('Crop_Eq')['rev_prior_yr'].to_dict()
    
    # Function to get a specific coefficient
    def pull_coef(df, cropname, varname):
        coef = df.loc[df['Crop_Eq'] == cropname, varname].values[0]
        return coef
    
    # Calculate predicted probabilities
    # pasture is the base category that fills the probability toward 1
    data['N.Mixed_pasture'] = 1
    for crop in crops:
        data[f'N.{crop}'] = np.exp(
            coef_rev_prior_yr * rpy_vector[crop] +
            pull_coef(logit_coefs, crop, 'intercept') +
            pull_coef(logit_coefs, crop, 'acres') * data['acres'] +
            pull_coef(logit_coefs, crop, 'elev_mean') * data['elev_mean'] +
            pull_coef(logit_coefs, crop, 'slope_perc_mean') * data['slope_perc_mean'] +
            pull_coef(logit_coefs, crop, 'pod') * data['pod'] +
            pull_coef(logit_coefs, crop, 'wy_dc') * data['wy_dc'] +
            pull_coef(logit_coefs, crop, 'i_pod_wy') * (data['pod'] * data['wy_dc']) +
            pull_coef(logit_coefs, crop, 'dtwfa') * data['dtwfa'] +
            pull_coef(logit_coefs, crop, 'dtwfa2') * (data['dtwfa']**2)
        )
    
    # need to make sure crop selection follows the same column order used for cumulative probability
    crops_N = data.loc[:,data.columns.str.contains('N.')].columns.str.replace('N.','').values
    # crops_N = logit_coefs['Crop_Eq'].tolist()

    # Calculate Denominator
    data['Denominator'] = 1 + data[[f'N.{crop}' for crop in crops_N]].sum(axis=1)
        
    # Calculate predicted probabilities
    for crop in crops_N:
        data[f'PP.{crop}'] = data[f'N.{crop}'] / data['Denominator']

    return data

# %%

# simplified version where predict_crops is pulling from the same base script
def predict_crops_alt(data):

    crops_N = data.loc[:,data.columns.str.contains('N.')].columns.str.replace('N.','').values
    # Find the maximum probability
    data['maxprob'] = data[[f'PP.{crop}' for crop in crops_N]].max(axis=1)
    
    # Assign Crop Choice based on maximum probability
    conditions = [data[f'PP.{crop}'] == data['maxprob'] for crop in crops_N]
    data['Crop_Choice'] = np.select(conditions, crops_N, default='Unclassified_fallow')
    # a suggestion by the peer review group was that a farmer may not always
    # select the highest probability, thus a random number should be used to select
    # from the distribution to determine the crop used
    # wait for Yusuke to run test with validation before updating
    
    # MODFLOW/irrigation choice info: parcel_id and Crop_Choice
    # not going to extract probabilities because that should be reviewed/processed elsewhere
    return data[['parcel_id', 'Crop_Choice']]


# simplified version where predict_crops is pulling from the same base script
def predict_crops_alt_rand(data):

    crops_N = data.loc[:,data.columns.str.contains('N.')].columns.str.replace('N.','').values
    
    # calculate cumulative probability
    data_PP_cum = data.loc[:,data.columns.str.contains('PP.')].cumsum(axis=1)
    # create random generator with seed for reproducibility
    rng = np.random.default_rng(seed=42)
    # get random numbers from 0 to 1
    rand_arr = rng.random(size=len(data_PP_cum))
    # identify the bracket which the random number falls in
    # rand_loc = data_PP_cum.values > np.repeat(np.reshape(rand_arr, (-1,1)), data_PP_cum.shape[1], axis=1)
    # subtract the random number from the cumulative probability
    rand_loc = data_PP_cum.values - np.repeat(np.reshape(rand_arr, (-1,1)), data_PP_cum.shape[1], axis=1)
    # where negative the random number exceeded the category
    rand_loc = np.ma.masked_where(rand_loc<0, rand_loc)
    # then select the first positive difference
    crop_rand = np.argmin(rand_loc, axis=1)
    # need to make sure crop selection follows the same column order used for cumulative probability
    # crop_select = data_PP_cum.columns.str.replace('PP.','').values
    # then select the crops based on the crop list
    data['Crop_Choice'] = np.array(crops_N)[crop_rand]
    return data[['parcel_id', 'Crop_Choice']]


# %%

def predict_crops_rand(data, rev_prior_yr_df, logit_coefs, return_prob=False):
    ''' Predict crop choice based on input variables
    data: dataframe with variables as columns and rows for each parcel
        'parcel_id', 'year', 'crop_cat', 'acres', 'elev_mean',
       'slope_perc_mean', 'dtwfa', 'dtwsp', 'pod'
    rev_prior_yr_df: dataframe of revenue expected for each crop
    logit_coefs: dataframe with the logit coefficients for each variable in data for each crop
    '''
    # Get a vector of crop equations
    crops = logit_coefs['Crop_Eq'].tolist()

    # Add water year critical or dry indicator
    data['dtwfa2'] = data['dtwfa'] * data['dtwfa']
    data['i_pod_wy'] = data['pod'] * data['wy_dc']
    
    # Get the coefficient on prior year revenue - same across crops
    coef_rev_prior_yr = logit_coefs.loc[0, 'rev_prior_yr']
    
    # Revenue for the previous year by crop
    rpy_vector = rev_prior_yr_df.set_index('Crop_Eq')['rev_prior_yr'].to_dict()
    

    
    # Calculate predicted probabilities
    # pasture is the base category that fills the probability toward 1
    data['N.Mixed_pasture'] = 1
    for crop in crops:
        data[f'N.{crop}'] = np.exp(
            coef_rev_prior_yr * rpy_vector[crop] +
            pull_coef(logit_coefs, crop, 'intercept') +
            pull_coef(logit_coefs, crop, 'acres') * data['acres'] +
            pull_coef(logit_coefs, crop, 'elev_mean') * data['elev_mean'] +
            pull_coef(logit_coefs, crop, 'slope_perc_mean') * data['slope_perc_mean'] +
            pull_coef(logit_coefs, crop, 'pod') * data['pod'] +
            pull_coef(logit_coefs, crop, 'wy_dc') * data['wy_dc'] +
            pull_coef(logit_coefs, crop, 'i_pod_wy') * (data['pod'] * data['wy_dc']) +
            pull_coef(logit_coefs, crop, 'dtwfa') * data['dtwfa'] +
            pull_coef(logit_coefs, crop, 'dtwfa2') * (data['dtwfa']**2) +
            # # new coefficients coef of last crop multiplied by 0 or 1 of last crop
            pull_coef(logit_coefs, crop, 'Vineyards') * data['Vineyards'] +
            pull_coef(logit_coefs, crop, 'Miscellaneous_grain_and_hay') * data['Miscellaneous_grain_and_hay'] +
            pull_coef(logit_coefs, crop, 'Corn__sorghum_or_Sudan') * data['Corn__sorghum_or_Sudan'] +
            pull_coef(logit_coefs, crop, 'Alfalfa_and_alfalfa_mixtures') * data['Alfalfa_and_alfalfa_mixtures']  +
            pull_coef(logit_coefs, crop, 'Unclassified_fallow') * data['Unclassified_fallow'] +
            # should there also be one for mixed pasture???
            pull_coef(logit_coefs, crop, 'Other') * data['Other']
            
        )
    
    # has all crops listed whiel logit_coefs were missing for pasture
    # crops_N = rev_prior_yr_df.Crop_Eq.to_list()
    # need to make sure crop selection follows the same column order used for cumulative probability
    crops_N = data.loc[:,data.columns.str.contains('N.')].columns.str.replace('N.','').values
    # crops_N = logit_coefs['Crop_Eq'].tolist()
    # crops_N = logit_coefs['Crop_Eq'].tolist()
    # Calculate Denominator
    # originally developed where crops was not including pasture
    # data['Denominator'] = 1 + data[[f'N.{crop}' for crop in crops_N]].sum(axis=1) # old version
    # sum should include 1 + sum of crops not including pivot crop (pasture)
    data['Denominator'] = data[[f'N.{crop}' for crop in crops_N]].sum(axis=1) # new version

    # Calculate predicted probabilities
    for crop in crops_N:
        data[f'PP.{crop}'] = data[f'N.{crop}'] / data['Denominator']

    # rather than assuming pasture has N=1 it might be more logical
    # to calculate it as the difference of sum of PP and 1
    if return_prob:
        return data
    
    # a suggestion by the peer review group was that a farmer may not always
    # select the highest probability, thus a random number should be used to select
    # from the distribution to determine the crop used
    # wait for Yusuke to run test with validation before updating

    # calculate cumulative probability
    data_PP_cum = data.loc[:,data.columns.str.contains('PP.')].cumsum(axis=1)
    # create random generator with seed for reproducibility
    rng = np.random.default_rng(seed=42)
    # get random numbers from 0 to 1
    rand_arr = rng.random(size=len(data_PP_cum))
    # identify the bracket which the random number falls in
    # rand_loc = data_PP_cum.values > np.repeat(np.reshape(rand_arr, (-1,1)), data_PP_cum.shape[1], axis=1)
    # subtract the random number from the cumulative probability
    rand_loc = data_PP_cum.values - np.repeat(np.reshape(rand_arr, (-1,1)), data_PP_cum.shape[1], axis=1)
    # where negative the random number exceeded the category
    rand_loc = np.ma.masked_where(rand_loc<0, rand_loc)
    # then select the first positive difference
    crop_rand = np.argmin(rand_loc, axis=1)
    # need to make sure crop selection follows the same column order used for cumulative probability
    # crop_select = data_PP_cum.columns.str.replace('PP.','').values
    # then select the crops based on the crop list
    data['Crop_Choice'] = np.array(crops_N)[crop_rand]
    
    # MODFLOW/irrigation choice info: parcel_id and Crop_Choice
    # not going to extract probabilities because that should be reviewed/processed elsewhere
    return data[['parcel_id', 'Crop_Choice']]


