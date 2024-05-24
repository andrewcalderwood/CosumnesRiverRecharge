
import h5py



def read_crop_arr_h5(crop, h5_fn):
    # convert arrays of annual rates to hdf5 files individually
    with h5py.File(h5_fn, "r") as f:
        grp = f['array'] 
        # print(grp.keys())
        dset = grp[crop]
        arr = dset[:]
    return(arr)
        
    