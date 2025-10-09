
rem set m_nam=input_write_2014_2020
rem set m_nam=input_write_2014_2022_R20
rem set m_nam=input_write_2014_2022_R3
rem set m_nam=input_write_2014_2022_R4

set m_nam=input_write_2014_2022_R200

set input_name=static_input_data_no_p_o.xlsx
rem set input_name=static_input_data.xlsx


rem call activate geo_env
call activate geo_env2
python 03b_summarize_output.py %2 %m_nam% %input_name%
PAUSE