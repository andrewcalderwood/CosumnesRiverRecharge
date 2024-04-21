# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:52:58 2023

@author: cameron.speir
Updated to function 1/23/2024
"""


import pandas as pd
import numpy as np
import os


def predict_crops(data, rev_prior_yr_df, logit_coefs):
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
    
    # Calculate Denominator
    data['Denominator'] = 1 + data[[f'N.{crop}' for crop in crops]].sum(axis=1)
    
    f'N.{crop}'
    
    # Calculate predicted probabilities
    for crop in crops:
        data[f'PP.{crop}'] = data[f'N.{crop}'] / data['Denominator']
    
    # Find the maximum probability
    data['maxprob'] = data[[f'PP.{crop}' for crop in crops]].max(axis=1)
    
    # Assign Crop Choice based on maximum probability
    conditions = [data[f'PP.{crop}'] == data['maxprob'] for crop in crops]
    data['Crop_Choice'] = np.select(conditions, crops)
    
    # MODFLOW/irrigation choice info: parcel_id and Crop_Choice
    # not going to extract probabilities because that should be reviewed/processed elsewhere
    return data[['parcel_id', 'Crop_Choice']]


