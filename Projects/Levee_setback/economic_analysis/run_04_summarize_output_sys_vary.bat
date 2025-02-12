set m_nam=input_write_2014_2020
rem set m_nam=input_write_2014_2022_R20
rem set m_nam=input_write_2014_2022_R3
rem set m_nam=input_write_2014_2022_R4
rem call activate geo_env
call activate geo_env2
python 03b_summarize_output.py %1 %m_nam%
PAUSE