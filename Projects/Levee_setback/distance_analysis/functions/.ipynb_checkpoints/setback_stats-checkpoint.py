'''
setback_stats to do simple statistic plotting
'''
import pandas as pd


#  reference numbers for the paper
# sample 600 m for large and long (type 2)
def compare_hf_to_all(setback, rch_xs_all, rch_xs_hf):
    all_s_med = rch_xs_all.loc[setback].iloc[1,:-1].median()
    hf_s_med = rch_xs_hf.loc[setback].iloc[1,:-1].median()
    print('For the '+str(setback)+'m setback')
    print('HCPs account for %.1f %%' %(100*hf_s_med/all_s_med), end=',')
    print('For the median recharge of %.1f MCM' %all_s_med)
    return(None)

def stats_lines(df, ax):
    """ Dataframe with realization as columns and setbacks as index"""
    # plot quantiles on the line plot
    quart = df.quantile([.25,.75], axis=1)
    quart.transpose().plot(color='tab:blue', ax=ax, legend=False)
    median = df.quantile([.5], axis=1)
    median.transpose().plot(color='tab:green', ax=ax, legend=False)
    # calculate whiskers
    iqr = quart.loc[0.75]-quart.loc[0.25]
    # 1.5 x the whole interquartile range
    whisker = pd.DataFrame(quart.loc[0.75] + iqr*1.5)
    whisker[1] = quart.loc[0.25] - iqr*1.5
    # where whisker is greater than max or min set as max or min
    whisker.loc[whisker[0]>df.max(axis=1), 0] = df.max(axis=1)[whisker[0]>df.max(axis=1)]
    whisker.loc[whisker[1]<df.min(axis=1), 1] = df.min(axis=1)[whisker[1]<df.min(axis=1)]
    whisker.plot(color='black', ax=ax, legend=False)
