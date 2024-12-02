rem set m_nam=input_write_2000_2022
rem m_nam=input_write_2000_2022_R20
set m_nam=input_write_2000_2022_R3
call activate geo_env
python 03b_summarize_output.py %1 %m_nam%
PAUSE