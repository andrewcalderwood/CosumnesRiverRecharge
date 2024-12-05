set m_nam=input_write_2014_2022
rem set m_nam=input_write_2014_2022_R20
rem set m_nam=input_write_2014_2022_R3
rem set m_nam=input_write_2014_2022_R4
call activate geo_env
python 03b_summarize_output.py %1 %m_nam%
PAUSE