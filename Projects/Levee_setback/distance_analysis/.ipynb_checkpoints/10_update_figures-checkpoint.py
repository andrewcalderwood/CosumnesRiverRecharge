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
from os.path import basename, dirname, join, expanduser
import glob

import pandas as pd

# %%
usr_dir = expanduser('~')
doc_dir = join(usr_dir, 'Documents')
    
git_dir = join(doc_dir,'GitHub','CosumnesRiverRecharge')
# dir of all gwfm data
gwfm_dir = join(usr_dir,'Box','research_cosumnes','GWFlowModel')

# %%
proj_dir = join(gwfm_dir, 'Levee_setback', 'levee_seback_distance_analysis')

writing_dir = join(usr_dir, 'Box', 'Classwork', 'HYD298_Writing in science')
paper_dir = join(writing_dir, 'Levee setback spatial','frontiers_in_submission','figures')

# %%
# create a folder to save the converted images
os.makedirs(join(paper_dir, 'tifs'), exist_ok=True)

# identify PNGs to convert
pngs = pd.Series(os.listdir(paper_dir))
pngs= pngs[pngs.str.contains('.png')]
# pngs = glob.glob(join(paper_dir, '*.png'))

# %%
from PIL import Image


# %%
pngs.iloc[0]

# %%
for n in range(len(pngs)):
    # Open the PNG image
    image = Image.open(join(paper_dir, pngs.iloc[n]))
    
    # Convert the image to TIFF
    image.save(join(paper_dir, 'tifs', 
                    pngs.iloc[n].split('.')[0]+'.tif'), "TIFF", dpi=(300, 300))
