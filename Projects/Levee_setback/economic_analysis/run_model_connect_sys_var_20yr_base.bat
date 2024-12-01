rem conda init activates a conda shell which is what conda.bat does
set m_nam=input_write_2000_2022
set scenario_name=none
echo %m_nam%
rem CALL conda.bat activate geo_env
call activate geo_env
python 03_model_connect.py %1 %m_nam% %scenario_name%
PAUSE