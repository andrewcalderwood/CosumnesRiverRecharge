rem conda init activates a conda shell which is what conda.bat does

rem R3:6x vineyard MAR, R4:90/20 floodplain MAR, R20:pumping constraint

rem set m_nam=input_write_2014_2022
rem set m_nam=input_write_2014_2022_R20
rem set m_nam=input_write_2014_2022_R3
rem m_nam=input_write_2014_2022_R4


rem version with no p_o
rem set m_nam=input_write_2014_2022_R204
set m_nam=input_write_2014_2022_R200
rem set m_nam=input_write_2014_2022_R203

echo %m_nam%

set input_name=static_model_inputs_no_p_o.xlsx
rem set input_name=static_model_inputs.xlsx

echo %input_name%

rem dont' use CALL conda.bat activate geo_env
rem laptop
rem call activate geo_env
rem local desktop
call activate geo_env
python 03_model_connect.py %2 %m_nam% %input_name%
rem python 03_model_connect.py %1 %m_nam%
PAUSE