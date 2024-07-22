
import pandas as pd
import numpy as np

class SoilDataProcessor:
    def __init__(self, soil_ag_all, CN, crop_in):
        self.soil_ag_in = soil_ag_all
        self.CN = CN
        self.crop_in = crop_in

    def load_and_process_soil_data(self, crop):
        self.soil_ag_all = pd.merge(self.soil_ag_in, self.CN)
        self.soil_ag_all = self.soil_ag_all.merge(self.crop_in, left_on='UniqueID', right_on='parcel_id')
        self.soil_ag_all = self.soil_ag_all[~self.soil_ag_all.name.isna()]
        self.soil_crop = self.soil_ag_all[self.soil_ag_all.name == crop].copy()
        return self

class SoilData:
    def __init__(self, soil_crop):
        self.soil_ag = soil_crop
        
    # what if this code was moved to the pre-processing
    # then the cleaned soil data could be avialable up front and for post-processing
    def prep_soil_data(self, etc_arr, var_crops):
        """
        Given a dataframe of soil properties create object variables for reference
        """
        # dictionary to be filled with soil data for Soil Water Budget
        self.nfield = len(self.soil_ag)
    
        # # when soil_K_low is missing using a substitute of Ksat/10
        self.Ks = np.where(self.soil_ag.Ksat_Low==0, self.soil_ag.Ksat/10, self.soil_ag.Ksat_Low)
        self.por = self.soil_ag.Porosity.values/100
        self.eps = self.soil_ag.EPS.values
        self.CN = self.soil_ag.CN.values
        self.psdi =  self.soil_ag.PSDI.values
        # parameter for Mualem, van Genuchten
        self.m = self.psdi/(self.psdi+1)
        self.wc_f =  self.soil_ag.w3rdbar.values/100 #field capacity
        self.wc_wp =  self.soil_ag.w15bar.values/100 #wilting point 
    
        # calculate the minimum of soil/crop rooting depth
        zr = (var_crops.zr/12)*0.3048
        self.depth = np.min([self.soil_ag.SoilDepth.values, np.repeat(zr,len(self.soil_ag))], axis=0)
        # self.depth = soil_ag.SoilDepth.values
        
        # Calculate total available water in the root zone
        self.taw = (self.wc_f - self.wc_wp)*self.depth
    
        # for runoff, convert CN from fraction to CN
        self.Smax = (1000/self.CN) - 10
        
        # p_table22: Soil water depletion fraction for no stress, requires ETc in mm/day
        self.P = var_crops['p_table22'] + 0.04*((5-(etc_arr*1000))) # Calculate adjusted daily soil water depletion fraction for no stress
        self.raw = self.taw*self.P # Calculate readily available water in the root zone
        return(self)


    # def run_ex(self):
    #     # Load and process soil data for the crop
    #     soil_crop = self.soil_processor.load_and_process_soil_data(self.crop)

    #     # Create and process soil_ag data
    #     etc_arr = np.zeros((len(season.dates)))
    #     soil_data = self.soil_processor.process_soil_ag_data(soil_crop, var_crops,





