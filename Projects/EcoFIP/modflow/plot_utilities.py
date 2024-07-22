"""
functions to support the baseline modflow model and concept output comparison
"""
import pandas as pd


def plt_wb(wb, wb0, plt_cols, plt_labels, ax, scale=1E-6):
    for n, var in enumerate(plt_cols):
        wb[var].multiply(scale).plot(ax=ax[n], 
                                     # label=label_baseline, 
                                     legend=False)
        wb0[var].multiply(scale).plot(ax=ax[n], 
                                      # label=label_restoration, 
                                      legend=False)

        ax[n].set_ylabel(plt_labels[n])

        ax[n].ticklabel_format(style='plain', axis='y')
        ax[n].set_xlabel(None)
#         ax[n].set_yscale('log')


# +

def sig_fill(ttest_all, ax):
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    for n in ttest_all.month.unique():
        t_df = ttest_all.loc[(ttest_all.month==n)]
        if t_df.sig.values[0]==1:
            t_min = t_df.month.min()-pd.DateOffset(days=15)
            t_max = t_df.month.max()+ pd.DateOffset(days=15) #+pd.DateOffset(months=1)
            ax.fill_between([t_min, t_max ], 
                             ylim[1], ylim[0], 
                            
                            color='blue', edgecolor='none', alpha=0.1) # hatch is too busy
# step/interpolate doesn't effect the box location
# could also use 'where=y > threshold' which would avoid the need for a loop


# -

def plt_wb_diff(wb_diff, plt_cols, plt_labels, ax, color, scale = 1E-3):
    for n, var in enumerate(plt_cols):
        wb_diff[var].multiply(scale).plot(ax=ax[n], label='Difference', legend=False, color=color)
        ax[n].set_ylabel(plt_labels[n])

        ax[n].ticklabel_format(style='plain', axis='y')
        ax[n].set_xlabel(None)

def get_lak_head(hdobj, lak_idx):
    """
    Return the spatially averaged head for the maximum head at the input locations (idx)
    hdobj: flopy head object
    idx: list of tuples as (layer, row, column)
    """
    # get heads under the lake
    lak_ts = hdobj.get_ts(lak_idx)
    lak_ts_df = pd.DataFrame(lak_ts, columns=['totim']+lak_idx)
    lak_ts_df = lak_ts_df.set_index('totim')
    lak_ts_df = lak_ts_df.melt(ignore_index=False)
    lak_ts_df[['k','i','j']] = lak_ts_df.variable.tolist()
    lak_ts_df = lak_ts_df.drop(columns='variable') # drop to speed up groupby
    lak_head = lak_ts_df.groupby(['totim','i','j']).max().groupby('totim').mean()
    return lak_head
