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

# %%
loadpth = 'C:/WRDAPP/GWFlowModel/Cosumnes/Regional/'
model_ws = loadpth+'historical_simple_geology_reconnection'

m = flopy.modflow.Modflow.load('MF.nam', model_ws=model_ws, 
                                exe_name='mf-owhm.exe', version='mfnwt')


# %%

# %% [markdown]
# # Load SWB output

# %%
def load_run_swb(crop, year):


    # %%
    def crop_arr_to_h5(arr, crop, h5_fn):
        # convert arrays of annual rates to hdf5 files individually
        with h5py.File(h5_fn, "a") as f:
            grp = f.require_group('array') # makes sure group exists
            grp.attrs['units'] = 'meters/day'
            grp.attrs['description'] = 'Rows represent the soil units and columns represent the days in the season'
            dset = grp.require_dataset(crop, arr.shape, dtype='f', compression="gzip", compression_opts=4)
            dset[:] = arr

# %% [markdown]
# # WEL update

# %%
fn = join(model_ws, 'field_SWB', "GW_applied_water_WY"+str(year)+".hdf5")
        with h5py.File(h5_fn, "r") as f:


# %%
# also need shapefile of pumping well locations for each parcel
parcel_wells = gpd.read_file(join(gwfm_dir, 'WEL_data', 'parcels_to_wells', 'parcels_to_wells.shp'))

# %%
# fields = pd.read_csv(join(uzf_dir, 'clean_soil_data', 'fields_output_reference.csv'))
# # extract applied water estimates where we have a known irrigaton type
# AW_irr = AW[:, fields.irr_name != 'no irrig']
fields_irr = fields[fields.irr_name != 'no irrig']

# %%
# identify the well to the field area
fields_well = pd.merge(fields_irr[['geom_id','name','irr_name','field_area_m2']],
         lu_wells[['geom_id','name','irr_name','row','column', 'depth_m', 'dist_m']], 
          how='left')

fields_spd = fields_well[['depth_m','row','column','field_area_m2']].copy()
fields_spd[['row','column']] -=1
frow = fields_spd.row
fcol = fields_spd.column
fields_spd['layer'] = get_layer_from_elev(dem_data[frow,fcol] - fields_spd.depth_m*0.9, botm[:, frow,fcol], m.dis.nlay)


# %% [markdown]
# # RCH update

# %%
fn = join(model_ws, 'field_SWB', "percolation_WY"+str(year)+".hdf5")


# %%
