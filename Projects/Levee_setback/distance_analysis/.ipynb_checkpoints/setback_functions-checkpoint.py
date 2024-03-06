

## legend and function to add quartiles, whiskers and median to line plot
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

stats_elements = [
    Line2D([0], [0],color='grey',label='Individual\nRealization'),
#     Line2D([0], [0],color='black',label='5th/95th', linestyle='--'),
    Line2D([0], [0],color='black',label='1.5x Quartile\nRange'),
    Line2D([0], [0],color='tab:blue',label='25th/75th'),
    Line2D([0], [0],color='tab:green',label='Median'),
]


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
#     quant = df.quantile([.05,.95], axis=1)
#     quant.transpose().plot(color='black', ax=ax, legend=False, linestyle='--')
