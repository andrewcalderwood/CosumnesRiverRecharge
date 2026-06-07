# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# Goal: script to help aggregate/clean files for sharing.
# - the head file size may be too large so that they may need to be simplified further

# %%
# standard python utilities
import os
import sys
from os.path import basename, dirname, join, exists, expanduser
import glob
import pandas as pd
import numpy as np
import time


# %%
import zipfile

# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir, 'Documents')
# dir of all gwfm data
# gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'


# %%
git_dir = join(doc_dir,'GitHub')
def add_path(fxn_dir):
    """ Insert fxn directory into first position on path so local functions supercede the global"""
    if fxn_dir not in sys.path:
        sys.path.insert(0, fxn_dir)
# flopy github path - edited
add_path(join(git_dir,'flopy'))
import flopy


# %%

loadpth =  'D:/WRDAPP/GWFlowModel'

loadpth = join(loadpth, 'Cosumnes/Stream_seepage')


# %%
upscale = 4 
upscale_txt = 'upscale'+str(upscale)+'x_'
model_nam = 'oneto_denier_'+upscale_txt+'2014_2018'

base_model_ws = join(loadpth,model_nam)
all_model_ws = join(loadpth, 'parallel_'+model_nam)


# %%

for r in np.arange(0,100): #100
    folder = 'realization'+ str(r).zfill(3)
    # update model workspace so outputs to right directory
    model_ws = join(all_model_ws, folder)

# %%
with zipfile.ZipFile(join(all_model_ws, 'MF_heads.zip'), 'w') as zip:
    for r in np.arange(0,10): #100
        folder = 'realization'+ str(r).zfill(3)
        # update model workspace so outputs to right directory
        model_ws = join(all_model_ws, folder)
        zip.write(join(model_ws, 'MF.hds'))
