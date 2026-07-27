rem conda init activates a conda shell which is what conda.bat does

set input_name=static_model_inputs_no_p_o.xlsx
rem set input_name=static_model_inputs.xlsx

echo %input_name%

rem R3:6x vineyard MAR, R4:90/20 floodplain MAR, R20:pumping constraint

rem version with no p_o and 2019 variables
set m_nam=input_write_2014_2025_R300
rem set m_nam=input_write_2014_2025_R304
rem use 303 for test with exist_shared
set m_nam=input_write_2014_2025_R303

rem version no p_o and use crop/SWB from R200
rem set m_nam=input_write_2014_2025_R404

echo %m_nam%

rem for variable pricing we use false
rem for model versions 300 we use 2019 values
set year_load_var_in=2019

echo %year_load_var_in%

rem turn on or off optimization to create rep_crop_soilbudget
set create_rep_swb=new_local
rem set create_rep_swb=existing_local
rem set create_rep_swb=existing_local_and_local_crop
set create_rep_swb=existing_shared

echo %create_rep_swb%

call activate geo_env
python 03_model_connect.py %4 %m_nam% %input_name% %year_load_var_in% %create_rep_swb%

PAUSE