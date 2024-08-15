from importlib import reload

import mf_model_class_test
reload(mf_model_class_test)
from mf_model_class_test import CosumnesWS, CosumnesModel

# initiate workspace with paths
py_ws = CosumnesWS()
# add relevant folder to path for loading modules
py_ws.setup_environment()

cosumnes = CosumnesModel(py_ws)

# may want to convert the model domain loading to be more interactive
cosumnes#.load_model_domain()
cosumnes.parent.gwfm_dir

# Example usage:
usr_dir = '~'
gwfm_dir = '/path/to/gwfm_dir'

model = CosumnesModel(usr_dir, gwfm_dir)
model.setup_environment()
xul, yul = model.load_model_domain()
model.define_model_attributes(ss_bool=False, strt_date='2000-10-01', end_date='2022-09-30')
model.configure_paths(scenario='reconnection')

gw_model = model.create_groundwater_model()
gw_model.configure_discretization(nrow=100, ncol=230, nlay=148, delr=200, delc=200, rotation=52.9, xul=xul, yul=yul, proj4_str='+proj=utm +zone=10 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')
dem_data = gw_model.load_elevation_data(gwfm_dir)
