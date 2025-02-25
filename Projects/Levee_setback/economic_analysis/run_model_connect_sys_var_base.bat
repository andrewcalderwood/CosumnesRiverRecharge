rem conda init activates a conda shell which is what conda.bat does

rem set m_nam=input_write_2014_2022
rem set m_nam=input_write_2014_2022_R20
rem set m_nam=input_write_2014_2022_R3
rem m_nam=input_write_2014_2022_R4


set m_nam=input_write_2014_2020
rem m_nam=input_write_2014_2020_R3
rem m_nam=input_write_2014_2020_R4


echo %m_nam%
rem CALL conda.bat activate geo_env
rem call activate geo_env
call activate geo_env2
python 03_model_connect.py %1 %m_nam% 
PAUSE