
import h5py
import numpy as np



def read_crop_arr_h5(crop, h5_fn):
    # convert arrays of annual rates to hdf5 files individually
    with h5py.File(h5_fn, "r") as f:
        grp = f['array'] 
        # print(grp.keys())
        dset = grp[crop]
        arr = dset[:]
    return(arr)


# +
def init_h5(h5_fn):
    """ Initiate hdf5 files for the given year before appending data for each crop"""
    with h5py.File(h5_fn, "w") as f:
        grp = f.require_group('array') # makes sure group exists
        # grp.attrs['units'] = 'meters/day'
        grp.attrs['description'] = 'Rows represent the soil units and columns represent the days in the season'

# these hdf5 files are written at the start so they can be appended to for each year and crop
def crop_arr_to_h5(arr, crop, h5_fn, units='meters/day'):
    # convert arrays of annual rates to hdf5 files individually
    with h5py.File(h5_fn, "a") as f:
        grp = f.require_group('array') # makes sure group exists
        grp.attrs['units'] = units
        # grp.attrs['description'] = 'Each layer of the array is a day in the water year'
        dset = grp.require_dataset(crop, arr.shape, dtype='f', compression="gzip", compression_opts=4)
        dset[:] = arr
