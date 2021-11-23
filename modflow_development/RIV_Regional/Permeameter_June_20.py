# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 16:17:16 2020

@author: arodri44
"""

#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# This code was created for Alisha to assist with the permeameter test with the MODFLOW model 
# to run various upscaling scenarios of TPROGS data


# In[ ]:


# Import needed packages

import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import flopy
import shapefile as sf
from flopy.utils.gridgen import Gridgen
from flopy.utils.reference import SpatialReference
import pandas as pd
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable
from flopy.utils.reference import SpatialReference
import flopy.utils.binaryfile as bf
from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'inline')
import scipy.io as spio


# Print the location and verion of flopy, numpy, matplotlib
print(sys.version)
print('numpy version: {}'.format(np.__version__))
print('matplotlib version: {}'.format(mpl.__version__))
print('flopy version: {}'.format(flopy.__version__))


# In[ ]:


modelname='Permeameter_June_2020'
#model_ws= '../Model/'
mf = flopy.modflow.Modflow(modelname, exe_name=r"C:\wrdapp\ucode_2014_1.004_and_more\ucode_2014_1.004\bin\mf2005.exe", version="mf2005")#model_ws=model_ws)


# In[ ]:


# Define domain and grid variables, DIS package
ztop = 0. # model top, msl
zbot = -100. # model bottom, msl
nlay = 1 # number of layers
nrow=60 # number of rows
ncol=105 # number of columns
delr=200 
delc=200 
laytyp = 0
Lx = delr * ncol
Ly = delc * nrow
delv = (ztop - zbot) / nlay
botm = np.linspace(ztop, zbot, nlay+1)
iter_mi=100
hclose =1E-2
rclose =1E-2
layavg =2
itmuni = 4. # indiates time unit of model data, 1 = seconds
# Define variables for BAS package
ibound = np.ones((nlay, nrow, ncol), dtype=np.int32) # defines active area of model
ibound[0, :, :] = -1 # constant head
ibound[-1, :, :] = -1 # constant head
strt = np.ones((nlay, nrow, ncol), dtype=np.float32)
strt[0, :, :] = 100. # defining constant head in starting array
strt[-1, :, :] = 90. # defining constant head in starting array
#strt = -5. * np.ones((nlay, nrow, ncol), dtype=np.float32)
#bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt)

 # Time step parameters
nper = 20
perlen = 100
nstp = 100
steady = False

# # Rotation and projection parameters
#xul=661602.879 # x, upper left hand corner of model grid
#yul=4272487.620 # y, upper left hand corner of model grid
#rotation=-47 # rotation
#proj4_str='EPSG:26910' # projection system


# In[ ]:


# Create the discretization object
dis = flopy.modflow.ModflowDis(mf, nlay, nrow, ncol, delr=delr, delc=delc,
                               top=ztop, botm=botm[1:],
                               nper=nper, perlen=perlen, nstp=nstp, steady=steady),
                               #xul=xul, yul=yul,
                               #rotation=rotation, proj4_str=proj4_str)
bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt)
pcg = flopy.modflow.ModflowPcg(mf,iter_mi, hclose=hclose, rclose=rclose)

# In[ ]:


# Instead of what is above, you need figure out what format your upscaled hk and vka, etc. 
# values are in and then import them by layer. For example:

vka_layers = np.ones((nlay, nrow, ncol), dtype=np.float32) # vertical K
hk_layers = np.ones((nlay, nrow, ncol), dtype=np.float32) # horizontal K
ss_layers = np.ones((nlay, nrow, ncol), dtype=np.float32) # Ss
sy_layers = np.ones((nlay, nrow, ncol), dtype=np.float32) # Sy

mat_vka = spio.loadmat('C:\Python\Kz_Upscale_4.mat')['Upscale_4']
mat_hk = spio.loadmat('C:\Python\K_xy_Upscale_4.mat')['Upscale_4']
mat_ss = spio.loadmat('C:\Python\Ss_Upscale_4.mat')['Upscale_4']
mat_sy = spio.loadmat('C:\Python\Sy_Upscale_4.mat')['Upscale_4']

# Then loop through geologic parameters (HK, VKA, SY, SS) for each layer
for l in range(nlay):
    vka_layers[l,:,:] = mat_vka[:,:,l] 
    hk_layers[l,:,:] = mat_vka[:,:,l] 
    ss_layers[l,:,:] = mat_ss[:,:,l] 
    sy_layers[l,:,:] = mat_sy[:,:,l]

#layvka = 0 # 0 indicates VKA is vertical hydraulic conductivity, not the ratio of H:V hydraulic conductivity 

lpf = flopy.modflow.mflpf.ModflowLpf(mf, layavg = layavg, hk=hk_layers, vka=vka_layers, ss=ss_layers, sy=sy_layers, laytyp = laytyp)





# In[ ]:


# Define variables for Output Control package
stress_period_data = {(0, 0): ['print head', 'print budget', 'save head', 'save budget']}

# User below script with >1 stress period
# stress_period_data = {}
# for kper in range(nper):
#     for kstp in range(nstp[kper]):
#         stress_period_data[(kper, kstp)] = ['save head',
#                                             'save drawdown',
#                                             'save budget',
#                                             'print head',
#                                             'print budget']
oc = flopy.modflow.ModflowOc(mf, stress_period_data=stress_period_data,
                             compact=True)


# In[ ]:


# Write the model input files
mf.write_input()

# Run the model
success, mfoutput = mf.run_model(silent=True, pause=False, report=True)
if not success:
    raise Exception('MODFLOW did not terminate normally.')


# In[ ]:


#fig = plt.figure(figsize=(20,20))
#plt.subplot(1, 1, 1, aspect='equal')
#hds = bf.HeadFile(modelname + '.hds')
#head = hds.get_data(totim=1.0)
##levels = np.arange(1, 10, 1)
#extent = (delr / 2., Lx - delr / 2., Ly - delc / 2., delc / 2.)
#cs = plt.contour(head[0, :, :], extent=extent)#levels=levels, 
#plt.clabel(cs, inline=1, fontsize=14, fmt='%1.1f')
##plt.savefig('tutorial1a.png')
#
#
## In[ ]:
#
#
#mf.plot()
#mf.dis.top.plot()


# In[ ]: