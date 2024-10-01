# ---
# jupyter:
#   jupytext:
#     formats: py:percent
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
# Code archive


# %%
# id_cols = ['layer','row','column','segment','reach']
# sfrdf_all = sfrdf.join(sfrdf0.set_index(id_cols, append=True), on=['dt']+id_cols, rsuffix='0')
# sfrdf_all = pd.concat((sfrdf.assign(scenario='baseline'), sfrdf0.assign(scenario='restoration')))
# sfrdf_all = pd.concat((sfrdf0.assign(scenario='baseline'), sfrdf.assign(scenario='restoration')))

# %% 
## doesn't help the story much because seepage by itself is only important if it improves water quality
## temperature or flow availability
# plt_df = sfrdf_all.copy()
# # plt_df[~sfrdf_all.facies.isin(['Mud'])] = np.nan

# plt_df = plt_df.groupby(['WY','Total distance (m)','scenario']).mean(numeric_only=True)
# # start simple with just year by segment ,'month','facies'
# sns.relplot(plt_df, x='Total distance (m)',y='Qrech', 
#             col='WY', col_wrap=2, hue='scenario', 
#             kind='line'
# #             kind='scatter'
#            )