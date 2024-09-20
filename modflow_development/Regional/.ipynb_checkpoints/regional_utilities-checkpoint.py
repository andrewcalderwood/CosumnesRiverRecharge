"""
regional_utilities module. 
Different functions for modflow set up with general python functions
First iteration as a Module January 2024
Author: Andrew Calderwood
"""

import pandas as pd
import numpy as np


def clean_wb(model_ws, dt_ref):
    # load summary water budget
    wb = pd.read_csv(model_ws+'/flow_budget.txt', delimiter=r'\s+')

    wb['kstpkper'] = list(zip(wb.STP-1,wb.PER-1))
    wb = wb.merge(dt_ref, on='kstpkper').set_index('dt')

    # calculate change in storage
    wb['dSTORAGE'] = wb.STORAGE_OUT - wb.STORAGE_IN
    # calculate total gw flow, sum GHB, CHD
    wb['GW_OUT'] = wb.GHB_OUT + wb.CHD_OUT
    wb['GW_IN'] = wb.GHB_IN + wb.CHD_IN
    wb = wb.loc[:,~wb.columns.str.contains('GHB|CHD')]
    
    wb_cols = wb.columns[wb.columns.str.contains('_IN|_OUT')]
    wb_cols = wb_cols[~wb_cols.str.contains('STORAGE|IN_OUT')]
    wb_out_cols= wb_cols[wb_cols.str.contains('_OUT')]
    wb_in_cols = wb_cols[wb_cols.str.contains('_IN')]
    # only include columns with values used
    wb_out_cols = wb_out_cols[np.sum(wb[wb_out_cols]>0, axis=0).astype(bool)]
    wb_in_cols = wb_in_cols[np.sum(wb[wb_in_cols]>0, axis=0).astype(bool)]

    return(wb, wb_out_cols, wb_in_cols)



