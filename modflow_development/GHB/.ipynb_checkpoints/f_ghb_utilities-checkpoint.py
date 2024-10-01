''' ghb functions module
First iteraton: 2024-7-26
'''
import pandas as pd
import numpy as np

def prep_ghb_df(ghb_df, hk, top_botm, delr=200):
    """ Given rows and columns create GHB based on interpolated head levels"""
    nlay = top_botm.shape[0]-1
    # pull out head for rows and columns
    rows = ghb_df.row.values
    cols = ghb_df.column.values
    ghb_lay = ghb_df.layer.values
    ghb_hd = ghb_df.set_index(['row','column'])
    head = ghb_hd.loc[list(zip(rows, cols))].value.values
    # ghb_lay = get_layer_from_elev(head, botm[:,rows, cols], m.dis.nlay) #0-based

    df = pd.DataFrame(np.zeros((np.sum(nlay - ghb_lay),5)))
    df.columns = ['k','i','j','bhead','cond']
    # get all of the i, j,k indices to reduce math done in the for loop
    n=0
    nk = 0 
    for i, j in list(zip(rows,cols)):
        for k in np.arange(ghb_lay[nk], nlay):
            df.loc[n,'i'] = i
            df.loc[n,'j'] = j
            df.loc[n,'k'] = k
            n+=1
        # update layer sampling location
        nk +=1

    df[['k','i','j']] = df[['k','i','j']].astype(int)
#     hk = hk[df.k, df.i, df.j] # old hk with cell by cell values
    distance = ghb_hd.loc[list(zip(df.i, df.j))].ghb_dist.values
    cond = hk*(top_botm[df.k, df.i, df.j]-top_botm[df.k +1 , df.i, df.j])*delr/distance
    df.cond = cond
    df.bhead = ghb_hd.loc[list(zip(df.i, df.j))].value.values
    # drop cells where the head is below the deepest cell?
    return(df)