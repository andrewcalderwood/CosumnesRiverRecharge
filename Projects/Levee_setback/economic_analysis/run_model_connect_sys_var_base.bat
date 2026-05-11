rem conda init activates a conda shell which is what conda.bat does

rem R3:6x vineyard MAR, R4:90/20 floodplain MAR, R20:pumping constraint

rem version with no p_o
set m_nam=input_write_2014_2025_R200
set m_nam=input_write_2014_2025_R203
rem set m_nam=input_write_2014_2025_R204


echo %m_nam%

set input_name=static_model_inputs_no_p_o.xlsx
rem set input_name=static_model_inputs.xlsx

echo %input_name%

call activate geo_env
python 03_model_connect.py %2 %m_nam% %input_name%
rem python 03_model_connect.py %1 %m_nam%
PAUSE