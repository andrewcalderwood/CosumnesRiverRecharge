import numpy as np
import pandas as pd
from os.path import join


class HydrogeologyModel:
    def __init__(self, CosumnesModel, botm, K, Ss, Sy, porosity, tprogs_info, upscale, num_leveling_layers, drain_layer, nlay_tprogs, grid_p, masked_tprogs, params, deep_geology, nlay, nrow, ncol, model_ws):
        self.botm = botm
        self.K = K
        self.Ss = Ss
        self.Sy = Sy
        self.porosity = porosity
        self.tprogs_info = tprogs_info
        self.upscale = upscale
        self.num_leveling_layers = num_leveling_layers
        self.drain_layer = drain_layer
        self.nlay_tprogs = nlay_tprogs
        self.grid_p = grid_p
        self.masked_tprogs = masked_tprogs
        self.params = params
        self.deep_geology = deep_geology
        self.nlay = nlay
        self.nrow = nrow
        self.ncol = ncol
        # flopy has the object model passed but the code still references parent
        # The Child class inherits from the Parent class using the syntax class Child(Parent)
        # An alternative approach (and not the only one) is to pass the parent as an argument to the child to create a property that corresponds to the parent:
        # nrow, ncol, nlay, nper = self.parent.nrow_ncol_nlay_nper
        self.model_ws = model_ws
        
        self.hk = np.zeros(self.botm.shape)
        self.vka = np.zeros(self.botm.shape)
        self.sy = np.zeros(self.botm.shape)
        self.ss = np.zeros(self.botm.shape)
        
    def hmean(self, arr, axis):
        return np.mean(arr, axis=axis)  # Modify this with the appropriate harmonic mean calculation if necessary

    def get_tprogs_for_elev(self, var, top, bot, tprogs_info):
        return tc.get_tprogs_for_elev(var, top, bot, tprogs_info)  # Assuming this function exists and is imported

    def process_layers(self):
        top = np.copy(self.botm[0, :, :])
        bot1 = np.copy(self.botm[-3, :, :])

        K_c = self.get_tprogs_for_elev(self.K, top, bot1, self.tprogs_info)
        Ss_c = self.get_tprogs_for_elev(self.Ss, top, bot1, self.tprogs_info)
        Sy_c = self.get_tprogs_for_elev(self.Sy, top, bot1, self.tprogs_info)
        n_c = self.get_tprogs_for_elev(self.porosity, top, bot1, self.tprogs_info)

        tprog_strt_lay = self.num_leveling_layers + self.drain_layer

        for kt, k in enumerate(np.arange(tprog_strt_lay, self.nlay_tprogs + tprog_strt_lay)):
            self.hk[k, :] = np.mean(K_c[self.upscale * kt:self.upscale * (kt + 1)], axis=0)
            self.vka[k, :] = self.hmean(K_c[self.upscale * kt:self.upscale * (kt + 1)], axis=0)
            self.ss[k, :] = np.mean(Ss_c[self.upscale * kt:self.upscale * (kt + 1)], axis=0)
            self.sy[k, :] = np.mean(Sy_c[self.upscale * kt:self.upscale * (kt + 1)], axis=0)

        top = m.dis.top.array
        bot1 = m.dis.botm.array[self.drain_layer, :, :]

        for k in np.arange(0, tprog_strt_lay):
            self.hk[k, :, :] = np.mean(self.get_tprogs_for_elev(self.K, top, bot1, self.tprogs_info), axis=0)
            self.vka[k, :, :] = self.hmean(self.get_tprogs_for_elev(self.K, top, bot1, self.tprogs_info), axis=0)
            self.sy[k, :, :] = np.mean(self.get_tprogs_for_elev(self.Sy, top, bot1, self.tprogs_info), axis=0)
            self.ss[k, :, :] = np.mean(self.get_tprogs_for_elev(self.Ss, top, bot1, self.tprogs_info), axis=0)

    def adjust_vka_quants(self):
        rows, cols = self.grid_p.row.values - 1, self.grid_p.column.values - 1

        tprogs_vals = np.arange(1, 5)
        tprogs_hist = np.histogram(self.masked_tprogs, np.append([0], tprogs_vals + 0.1))[0]
        tprogs_hist = tprogs_hist / np.sum(tprogs_hist)

        tprogs_quants = 1 - np.append([0], np.cumsum(tprogs_hist) / np.sum(tprogs_hist))
        vka_quants = pd.DataFrame(tprogs_quants[1:], columns=['quant'], index=tprogs_vals)

        vka_quants['vka_min'] = np.quantile(self.vka, tprogs_quants[1:])
        vka_quants['vka_max'] = np.quantile(self.vka, tprogs_quants[:-1])
        vka_quants['facies'] = self.params.loc[tprogs_vals].Lithology.values

        for p in tprogs_vals:
            self.vka[(self.vka < vka_quants.loc[p, 'vka_max']) & (self.vka > vka_quants.loc[p, 'vka_min'])] /= self.params.vani[p]

        vka_quants.to_csv(join(self.model_ws, 'vka_quants.csv'))

    def set_laguna_and_mehrten_layers(self):
        self.hk[-2, :, :] = self.params.loc[5, 'K_m_d']
        self.vka[-2, :, :] = self.params.loc[5, 'K_m_d'] / self.params.loc[5, 'vani']
        self.sy[-2, :, :] = self.params.loc[5, 'Sy']
        self.ss[-2, :, :] = self.params.loc[5, 'Ss']

        self.hk[-1, :, :] = self.params.loc[6, 'K_m_d']
        self.vka[-1, :, :] = self.params.loc[6, 'K_m_d'] / self.params.loc[6, 'vani']
        self.sy[-1, :, :] = self.params.loc[6, 'Sy']
        self.ss[-1, :, :] = self.params.loc[6, 'Ss']

    def adjust_deep_geology(self):
        adj_lowK = pd.DataFrame(np.transpose(np.where(self.deep_geology > 0)), columns=['k', 'i', 'j'])
        adj_lowK = adj_lowK.groupby('k').quantile(0.15)['j'].astype(int)

        adj_lowK_arr = np.zeros((self.nlay, self.nrow, self.ncol))
        for k in adj_lowK.index:
            adj_lowK_arr[k, :, adj_lowK.loc[k]:] = 1

        adj_lowK_arr = adj_lowK_arr.astype(bool)

        param_d = 6
        self.hk[adj_lowK_arr] = self.params.loc[param_d, 'K_m_d']
        self.vka[adj_lowK_arr] = self.params.loc[param_d, 'K_m_d'] / self.params.loc[param_d, 'vani']
        self.sy[adj_lowK_arr] = self.params.loc[param_d, 'Sy']
        self.ss[adj_lowK_arr] = self.params.loc[param_d, 'Ss']

        if self.drain_layer == 1:
            self.hk[0, adj_lowK_arr[0]] = self.params.loc[5, 'K_m_d']
            self.vka[0, adj_lowK_arr[0]] = self.params.loc[5, 'K_m_d'] / self.params.loc[5, 'vani']

    def run_model(self):
        self.process_layers()
        self.adjust_vka_quants()
        self.set_laguna_and_mehrten_layers()
        self.adjust_deep_geology()

# Example usage:
# model = HydrogeologyModel(botm, K, Ss, Sy, porosity, tprogs_info, upscale, num_leveling_layers, drain_layer, nlay_tprogs, grid_p, masked_tprogs, params, deep_geology, nlay, nrow, ncol, model_ws)
# model.run_model()
