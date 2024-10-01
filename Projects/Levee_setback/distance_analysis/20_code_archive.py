## Code chunk archive


## old output loading


def load_output(ft_in, region):
    tic = time.time()
    T_in = int(10**flood_type.loc[ft_in,'log_no_d'])
    p_l_in = flood_type.loc[ft_in,'pk_loc']
    tp_in = int(p_l_in*T_in)
    rch_hf_all = np.zeros((100, len(setbacks),nrow,ncol))
    Q_all = np.zeros((100, T_in, len(setbacks),xs_levee_smooth.shape[1]+1))

    # filter out for only those realizations that successfully ran
    base_fn = join(data_dir, region, 'type'+str(ft_in))
    r_out = pd.Series(os.listdir(base_fn)).str.extract(r'(\d{3})')[0].unique().astype(int)
    # takes a 
    for t in r_out: # np.arange(0,100): #[0]:
        r_fn = join(base_fn,'r'+str(t).zfill(3)+'_')
        # saving all of the flow at all steps, setbacks is needed to post-process
        Q_in = np.loadtxt(r_fn+'flow.tsv')
        Q = np.reshape(Q_in, ((T_in, len(setbacks), xs_levee_smooth.shape[1]+1)))
        Q_all[t,:] = np.copy(Q)

        # for recharge we want to aggregate across time steps but look at differences across setbacks
        rch_in = np.loadtxt(r_fn+'recharge.tsv')
        rch_sum = np.reshape(rch_in, (len(setbacks), nrow, ncol))
        rch_hf_all[t] = np.copy(rch_sum)
    # convert to m3/day and will have the total recharged after summing individual days
    rch_hf_all = rch_hf_all*86400

    # convert to total AF from total m3
    # rch_hf_all = rch_hf_all/((0.3048**3)*43560)

    toc = time.time()
    print('Loading',region,'for flow type',str(ft_in), 'took %.2f minutes' %((toc-tic)/60))
    return(Q_all, rch_hf_all)


## old post-calculation of depth from flow

# rating curves for each segment and setback
xs_flow_all = pd.read_csv(join(chan_dir,'all_xs_50pt_rating_curves.csv'))

def depth_match(seg_flow, flow):
    """ Given a XS (nseg, setback) return the expected depth (m) given a flow (cms)"""
    # find flows above and below the input flow
    flow_diff = (seg_flow.flow_cms-flow)
    f_high = flow_diff[flow_diff>0].argsort().index[0]
    f_low = flow_diff[flow_diff<0].argsort().index[-1]
    match_d = seg_flow.loc[[f_low, f_high]].sort_values('flow_cms')
    # linearly interpolate to calculate exact depth
    flow_slope = (match_d.iloc[1].flow_cms-match_d.iloc[0].flow_cms)/(match_d.iloc[1].depth_m-match_d.iloc[0].depth_m)
    out_depth = match_d.iloc[0].depth_m + (flow-match_d.iloc[0].flow_cms)/flow_slope
    return(out_depth)
    
# # this method is slightly different than the depth from the output because here
# the depth from the XS is given vs the depth above the raster
# this method which back calculates the depth given the flow
# returns the maximum depth for each segment which is what fish would likely have access to??
# recreate_hdf5 = True
recreate_hdf5 = False

if recreate_hdf5:
    t0 = time.time()
    region = 'regional'
    for ft_in in [1,2,3]:
        f = h5py.File(join(data_dir,'hdf5', 'all_flow_'+region+'_type'+str(ft_in)+'.hdf5'), "r")
        Q_all = f['array']['all'][:]
        f.close()
        d_arr = np.zeros(Q_all.shape)
        # iterate over realizations
        for r in np.arange(Q_all.shape[0]):
            # days of flow
            for t in np.arange(0, Q_all.shape[1]):
                # setbacks
                for s, setback in enumerate(setbacks):
                    # segments
                    for nseg in np.arange(0,num_segs):
                        seg_flow = xs_flow_all[(xs_flow_all.nseg==nseg)&(xs_flow_all.setback==setback)]
                        if Q_all[r,t, s, nseg] >=seg_flow.flow_cms.min():
                            d_arr[r,t, s, nseg] = depth_match(seg_flow, flow=Q_all[r,t, s, nseg])
                        else:
                            d_arr[r,t, s, nseg] = 0
        # save to hdf5 file for each flow type
        arr_to_h5(d_arr, join(data_dir,'hdf5', 'peak_flow_xs_depth_'+region+'_type'+str(ft_in)+'.hdf5'))

    # output will be all depths
    t1 = time.time()
    # 6.7 seconds for 1 realization, ~10 min for 100 x 3 flood types is ~30 min
    # the shortest flow type probably only took 10 min, the others are like 5-10x longer
    # seems to take longer than 30 min because it took 2 hrs in total
    print((t1-t0)/60)
else:
    print('Reusing existing recharge and flow hdf5 files')