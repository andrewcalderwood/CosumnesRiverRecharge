rem conda init activates a conda shell which is what conda.bat does

set input_name=static_model_inputs_no_p_o.xlsx

echo %input_name%

rem for model versions 300 we use 2019 values
set year_load_var_in=2019
rem for variable pricing we use false
rem set year_load_var_in=False

echo %year_load_var_in%


call activate geo_env
python 01_rep_SWB_runs.py %2 %input_name% %year_load_var_in%
PAUSE