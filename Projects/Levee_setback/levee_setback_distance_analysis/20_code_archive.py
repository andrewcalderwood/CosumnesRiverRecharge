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