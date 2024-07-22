import os
import time
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.optimize import Bounds

# class created from chatgpt to test utility
# it seems like it would help with variable referencing but take a bit of re-working

class CropModel:
    def __init__(self, crop, year, base_model_ws, proj_dir, uzf_dir, dtw_df, soil_rep=False, run_opt=True):
        self.crop = crop
        self.year = year
        self.base_model_ws = base_model_ws
        self.proj_dir = proj_dir
        self.dtw_df = dtw_df
        self.soil_rep = soil_rep
        self.run_opt = run_opt
        
        self.var_gen, self.var_crops, self.var_yield, self.season, self.pred_dict, self.crop_dict = swb.load_var(crop)
        self.yield_start = swb.ymd2dt(year, self.season.month_start, self.season.day_start, self.season.start_adj)
        self.yield_end = swb.ymd2dt(year, self.season.month_end, self.season.day_end, self.season.end_adj)
        self.strt_date = self.yield_start.min()
        self.end_date = self.yield_end.max()
        self.dates = pd.date_range(self.strt_date, self.end_date, freq='D')
        self.nper = (self.end_date - self.strt_date).days + 1
        
        self.initialize_swb_folder()
        
        self.Kc, self.Kc_dates = swb.load_Kc(year)
        self.prepare_Kc_df()
        
        self.rain, self.ETo_df = swb.load_hyd(year, self.dates)
        self.ETc = self.ETo_df.values * self.Kc_df.Kc.values
        
        self.prepare_gen_dict()

        # the irrigation efficiency should be changed to a subclass
        self.avg_irr_eff = pd.read_csv(os.path.join(proj_dir, 'model_inputs', 'avg_irr_eff_by_crop.csv'), index_col=0)
        self.gen_dict['irr_eff_mult'] = 100 / self.avg_irr_eff.loc[self.crop_dict[self.crop]].Avg_eff
        
        self.soil_crop = swb.load_soil(self.pred_dict[self.crop], self.crop_dict[self.crop])
        if self.soil_rep:
            self.soil_crop = self.soil_crop.groupby('pod').mean(numeric_only=True).reset_index()
            self.soil_crop = pd.concat([self.soil_crop] * dtw_df.shape[1])
        self.nfield_crop = len(self.soil_crop)
        
        self.prepare_crop_dtw()
        
        self.gen = cost_variables(self.gen_dict)
        self.bounds = Bounds(lb=0)
        self.soil_df_out = pd.DataFrame()

    def initialize_swb_folder(self):
        os.makedirs(os.path.join(self.base_model_ws, 'field_SWB'), exist_ok=True)
        for var in ['profit', 'yield', 'percolation', 'GW_applied_water', 'SW_applied_water']:
            name = os.path.join(self.base_model_ws, 'field_SWB', f'{var}_WY{self.year}.hdf5')
            init_h5(name)
    
    def prepare_Kc_df(self):
        Kc_c = self.Kc[self.Kc.Crop == self.crop]
        Kc_dates_c = self.Kc_dates[self.Kc_dates.Crop == self.crop]
        
        self.Kc_df = pd.DataFrame()
        for c in Kc_dates_c.Cycle.unique():
            Kc_df_c = swb.get_Kc_dates(Kc_dates_c[Kc_dates_c.Cycle == c], Kc_c)
            self.Kc_df = pd.concat((self.Kc_df, Kc_df_c))
    
    def prepare_gen_dict(self):
        self.var_crops_dict = {v: self.var_crops[v].tolist() for v in self.var_crops.index.unique()}
        
        self.gen_dict = {**self.var_gen.to_dict(), **self.var_crops_dict}
        self.gen_dict['y_max'] = self.var_crops[['y_max']].values
        self.gen_dict['nper'] = self.nper
        self.gen_dict['yield_ind'] = np.append([0], (self.yield_end - self.strt_date).dt.days.values + 1)
        
        self.var_yield['dt'] = swb.ymd2dt(self.year, self.var_yield.month, self.var_yield.day, self.var_yield.year_adj)
        self.gen_dict['K_Y'] = self.var_yield.set_index('dt')['K_Y'].reindex(self.dates).ffill().values
        
        gap_irr = self.var_crops['gap_irr']
        self.gen_dict['n_irr'] = np.floor(len(self.dates) / gap_irr).astype(int) + 1
        self.gen_dict['irr_days'] = np.arange(0, (self.gen_dict['n_irr'] * gap_irr - 1), gap_irr).astype(int)
    
    def prepare_crop_dtw(self):
        if self.soil_rep:
            self.crop_dtw = self.dtw_df.copy()
            ntimes = int(self.nfield_crop / self.crop_dtw.shape[1])
            self.crop_dtw = pd.concat([self.crop_dtw] * ntimes, axis=1)
        else:
            self.crop_dtw = self.dtw_df.loc[:, self.soil_crop['UniqueID'].values]
        
        self.crop_dtw = self.crop_dtw.loc[self.dates].values
    
    def run_optimization(self):
        if not self.run_opt:
            return
        
        t0 = time.time()
        irr_all = np.zeros((self.nfield_crop, 2 * self.gen_dict['n_irr']))
        p_all = np.ones(self.nfield_crop)
        t_all = np.zeros(self.nfield_crop)
        
        for ns in range(self.nfield_crop):
            soil_ag = self.soil_crop.iloc[[ns]]
            nfield = soil_ag.shape[0]
            dtw_arr = self.crop_dtw[:, ns]
            water_source = choose_water_source(dtw_arr, self.gen, mix_fraction=1)
            soil_dict = prep_soil_dict(soil_ag, self.ETc, self.var_crops)
            soil = cost_variables(soil_dict)
            field_soil_df = pd.DataFrame(np.append([soil_dict[k][0] for k in soil_keys_keep], soil_ag.UniqueID.iloc[0])).transpose()
            self.soil_df_out = pd.concat((self.soil_df_out, field_soil_df), axis=0)
            
            sw_con = 100
            gw_con = 100
            if soil_ag.pod.iloc[0] == 'No Point of Diversion on Parcel' or water_source == 'gw':
                sw_con = 0
                n_irr_type = 1
            if water_source == 'sw':
                gw_con = 0
                n_irr_type = 1
            
            irr_lvl = np.ones(n_irr_type * self.gen_dict['n_irr']) * (2 / 12) * 0.3048
            if ns > 0:
                if water_source == 'gw':
                    irr_lvl[:] = irr_all[ns-1, self.gen_dict['n_irr']:]
                elif water_source == 'sw':
                    irr_lvl[:] = irr_all[ns-1, :self.gen_dict['n_irr']]
                else:
                    irr_lvl[:] = irr_all[ns-1]
            
            linear_constraint = mak_irr_con_adj(self.gen_dict['n_irr'], gw_con=gw_con, sw_con=sw_con)
            tol = 0.01
            
            while p_all[ns] > 0:
                out = minimize(run_swb, irr_lvl, args=(soil, self.gen, self.rain, self.ETc, dtw_arr, water_source),
                               method='trust-constr', constraints=[linear_constraint], bounds=self.bounds, tol=tol)
                if out.fun > 0:
                    tol /= 10
                    irr_lvl[:] = (2 / 12) * 0.3048
                if tol < 1E-5:
                    break
                
                if water_source == 'gw':
                    irr_all[ns, self.gen_dict['n_irr']:] = out.x
                elif water_source == 'sw':
                    irr_all[ns, :self.gen_dict['n_irr']] = out.x
                else:
                    irr_all[ns] = out.x
                
                p_all[ns] = out.fun
                t_all[ns] = out.execution_time
            
            print(f'Soil {ns} %.2f $ in %.2f min' % (-out.fun, out.execution_time / 60))
        
        t1 = time.time()
        print(f'Total time was %.2f min for {ns+1} parcels' % ((t1 - t0) / 60))
        return(irr_all, p_all)

# Usage
crop_model = CropModel(crop, year, base_model_ws, proj_dir, uzf_dir, dtw_df, soil_rep)
crop_model.run_optimization()
