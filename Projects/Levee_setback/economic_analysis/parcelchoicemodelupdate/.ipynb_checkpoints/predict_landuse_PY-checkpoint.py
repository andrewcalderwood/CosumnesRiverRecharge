# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:52:58 2023

@author: cameron.speir
"""

import pandas as pd
import numpy as np
import os

#

# +
# # Set working directory (not needed if input data are in folder)
# os.chdir('U:/SESYNC/Cosumnes_Ag')

# Read logit coefficients
logit_coefs = pd.read_csv('logit_coefs.csv')
# -

# Get a vector of crop equations
crops = logit_coefs['Crop_Eq'].tolist()

# the parcel data needs the dtwfa (avg dtw in feet for the parcel) and wy_dc (pulled from Sac wy type dataset and switched to dry boolean)
# missing WY type prediction? 
# Read parcel data
data = pd.read_csv("data_model/parcel_data_test.csv")

# Import prior year revenue data by crop
rev_prior_yr_df = pd.read_csv("data_model/rev_prior_yr.csv")

# Add water year critical or dry indicator
data['wy_dc'] = np.where(data['year'] == 2020, 1, 0) # should be pulled from Sac WY type data instead
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
data[['parcel_id', 'Crop_Choice']].to_csv('parcel_crop_choice.csv', index=False)


