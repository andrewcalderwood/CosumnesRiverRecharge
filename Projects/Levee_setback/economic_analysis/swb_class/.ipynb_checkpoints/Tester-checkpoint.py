# +
from os.path import join, exists, dirname, basename, expanduser
import pandas as pd

from importlib import reload

# +
usr_dir = expanduser('~')
doc_dir = join(usr_dir,'Documents')
# dir of all gwfm data
gwfm_dir = dirname(doc_dir)+'/Box/research_cosumnes/GWFlowModel'

uzf_dir = join(gwfm_dir,'UZF_data')

proj_dir = join(dirname(doc_dir),'Box','SESYNC_paper1')
data_dir = join(proj_dir, 'model_inputs')

# +
field_ids = 'parcel'
soil_path = join(uzf_dir,'clean_soil_data')
# soil data for each ag field
soil_ag_all = pd.read_csv(join(soil_path, 'soil_for_'+field_ids+'_fields.csv'), index_col=0)

# curve numbers
CN = pd.read_csv(join(soil_path, field_ids+'_field_CN.csv'),index_col = 0)


# +
crop_choice_dir = '../parcelchoicemodelupdate'
year=2020
data_out = pd.read_csv(join(crop_choice_dir, 'parcel_crop_choice_'+str(year)+'.csv'))
loadpth = 'C://WRDAPP/GWFlowModel/Cosumnes/Regional/'

base_model_ws = loadpth + 'crop_soilbudget'
crop_in = pd.read_csv(join(base_model_ws, 'field_SWB', 'crop_parcels_'+str(year)+'.csv'))

# -

import soil_class
reload(soil_class)

soil_processor = soil_class.SoilDataProcessor(soil_ag_all, CN, crop_in) 

# transform soil data into clean dataframe
soil_processor = soil_processor.load_and_process_soil_data('Vineyards')



etc_arr = np.zeros((nper))
for n in np.arange(nper):
    etc_arr[n] = ETc[n]

soil = soil_class.SoilData(soil_processor.soil_crop)
soil.prep_soil_data()


# +
class Tester:
    def __init__(self):
        self.base_name = 'base model'

class Tester_sub(Tester):
    def __init__(self, Tester):
        self.sub_name = 'sub 1'
        super(Tester_sub, self).__init__()


# -

basetest = Tester()

subtest = Tester_sub(basetest)
print(subtest.sub_name)
print(basetest.base_name)
print(subtest.base_name)
