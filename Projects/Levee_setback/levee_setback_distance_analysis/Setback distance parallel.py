## Setback distance parallizeable
import os
import flopy 
import pandas as pd  
import numpy as np  


## Set up directory referencing
# Package data
git_dir = os.path.dirname(os.path.dirname(os.getcwd()))
git_dir = os.getcwd()
while os.path.basename(git_dir) != 'CosumnesRiverRecharge':
    git_dir = os.path.dirname(git_dir)
gwfm_dir = '\\'.join(str.split(git_dir,'\\')[0:3])+ '/Box/research_cosumnes/GWFlowModel'
print(git_dir, gwfm_dir)

box_dir = gwfm_dir+'/Levee_setback/levee_setback_distance_analysis/'

# set base model path
ext_dir = 'F:/WRDAPP'
c_dir = 'C:/WRDAPP'

if os.path.exists(ext_dir):
    loadpth = ext_dir 
elif os.path.exists(c_dir):
    loadpth = c_dir 

loadpth = loadpth +'/GWFlowModel/Cosumnes/levee_setback/setback_distance_analysis/'
model_ws = loadpth+'Permeameter_for_velocity'

####################################################################################################
## Adjustments here ##
flow_percentile = 95
# tprogs_id = ''
tprogs_id = '_no_conditioning'

# data_dir = box_dir+'data_output/'
data_dir = box_dir+'no_conditioning/data_output/'


####################################################################################################
model_ws = model_ws +tprogs_id

# load the model
name = 'MF.nam'
m = flopy.modflow.Modflow.load(name, model_ws=model_ws, 
                                exe_name='mf2005', version='mf2005')


# Load grid and SFR grid for setback distances
## dem data for cropping above land surface
dem_data = np.loadtxt(gwfm_dir+'/DIS_data/dem_52_9_200m_linear.tsv')
grid_sfr = gpd.read_file(gwfm_dir+'/SFR_data/final_grid_sfr/grid_sfr.shp')
grid_p = gpd.read_file(gwfm_dir+'/DIS_data/grid/grid.shp')

# Calculate setback distances
buf_sfr = grid_sfr.copy()
setbacks = np.arange(0, 3400,200)
str_setbacks = np.zeros((len(setbacks),m.dis.nrow,m.dis.ncol))
# grid_sfr.plot()
for n in np.arange(0,len(setbacks)):
    buf_sfr.geometry = grid_sfr.buffer(setbacks[n])
    grid_sfr_buf = gpd.sjoin(grid_p,buf_sfr, how='right', lsuffix = 'grid', rsuffix = 'sfr',op='within')
    grid_sfr_buf = grid_sfr_buf.drop_duplicates('node_x')
    str_setbacks[n, grid_sfr_buf.row_x.values-1,grid_sfr_buf.column_x.values-1] = 1


def elev_to_tprogs_layers(elev, tprogs_top_elev, tprogs_bot_elev, num_lays):
    """
    function to get the tprogs layers based on the given elevation
    Example
    layer 0 is 80 meters, layer 1 is 79.5 meters, layer -1 is -80 meters
    """
    lay_thick = (tprogs_top_elev - tprogs_bot_elev)/num_lays
    elev_round = np.round((elev) * (1/lay_thick)) / (1/lay_thick) # dem rounded to the layer thickness
    elev_round[elev_round >= tprogs_top_elev] = tprogs_top_elev# any elevation above the top is set to the top
    # subtract the calculated row from top elev divided by layer thickness to get to index 0 at top and index 320 and bottom
    elev_indices = tprogs_top_elev/lay_thick - elev_round*(1/lay_thick) 
    return(elev_indices.astype(int))

# find tprogs layers at land surface for each row, column
tprogs_lay = elev_to_tprogs_layers(elev=dem_data,tprogs_top_elev=80, tprogs_bot_elev=-80, num_lays=320)
rows = np.where(np.ones(tprogs_lay.shape)==1)[0]
cols = np.where(np.ones(tprogs_lay.shape)==1)[1]

# will count total number of cells for each setback distance and for all 100 realizations
hf_tot = np.zeros(100,(len(setbacks)))

for r in np.arange(0,100):
	folder = '/realization'+ str(n).zfill(3)+'/'
	run_ws = model_ws+folder

	# iterable part
	cbb = flopy.utils.CellBudgetFile(run_ws+'/MF.cbc')
	# load velocity in z direction
	extcbb = flopy.utils.postprocessing.get_extended_budget(cbb)

	(qx, qy, qz) = flopy.utils.postprocessing.get_specific_discharge(vectors = extcbb, model=m)
	# convert flow to positive as it is all moving in the downward, -z direction
	qz *= -1

	# split cells into low and high conductivity, based on chosen flow percentile
	# shouldn't be less than 76th as only 24% of model is sand/gravel
	qz_plt = np.zeros(qz.shape)
	qz_plt[qz>=np.percentile(qz,[flow_percentile])] = 1

	# get high conductivity at ground surface
	qz_lay = np.zeros((100,230))
	qz_lay[rows,cols] = qz_plt[tprogs_lay[rows,cols],rows,cols]
	for n in np.arange(0,len(setbacks)):
	    hf = qz_lay[str_setbacks[n,:,:].astype('bool')]
	    hf_tot[r,n] = hf.sum()

hf_tot_df = pd.DataFrame(hf_tot, columns = setbacks)
# save counted high flow cells to a csv
hf_tot_df.to_csv(data_dir + 'surface_highflow_by_distance.csv')

