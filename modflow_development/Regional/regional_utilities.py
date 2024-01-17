"""
regional_utilities module. 
Different functions for modflow set up with general python functions
First iteration as a Module January 2024
Author: Andrew Calderwood
"""



def clean_sfr_df(model_ws, dt_ref, pd_sfr=None):
    sfrout = flopy.utils.SfrFile(join(model_ws, m.name+'.sfr.out'))
    sfrdf = sfrout.get_dataframe()
    sfrdf = sfrdf.join(dt_ref.set_index('kstpkper'), on='kstpkper').set_index('dt')
    # convert from sub-daily to daily using mean, lose kstpkper
    sfrdf = sfrdf.groupby(['segment','reach']).resample('D').mean(numeric_only=True)
    sfrdf = sfrdf.reset_index(['segment','reach'], drop=True)
    sfrdf[['row','column']]-=1 # convert to python
    sfrdf['month'] = sfrdf.index.month
    sfrdf['WY'] = sfrdf.index.year
    sfrdf.loc[sfrdf.month>=10, 'WY'] +=1
    # add column to track days with flow
    sfrdf['flowing'] = 1
    sfrdf.loc[sfrdf.Qout <= 0, 'flowing'] = 0
    if pd_sfr is not None:
    #     sfrdf = pd_sfr.join(sfrdf.set_index(['row','column']),on=['row','column'],how='inner',lsuffix='_all')
        sfrdf = sfrdf.join(pd_sfr ,on=['segment','reach'],how='inner',lsuffix='_all')

    # create different column for stream losing vs gaining seeapge
    sfrdf['Qrech'] = np.where(sfrdf.Qaquifer>0, sfrdf.Qaquifer,0)
    sfrdf['Qbase'] = np.where(sfrdf.Qaquifer<0, sfrdf.Qaquifer*-1,0 )
    if 'gradient' in sfrdf.columns:
        # booleans for plotting
        sfrdf['gaining'] = (sfrdf.gradient == 0)
        sfrdf['losing'] = (sfrdf.gradient >= 0)
        sfrdf['connected'] = (sfrdf.gradient < 1)
    return(sfrdf)

