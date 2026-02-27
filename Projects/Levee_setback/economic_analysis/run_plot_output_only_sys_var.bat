rem conda init activates a conda shell which is what conda.bat does


rem version with no p_o
set m_nam=input_write_2014_2022_R204
rem set m_nam=input_write_2014_2022_R200
rem set m_nam=input_write_2014_2022_R203

rem check to see if different start years have an impact
set m_nam=input_write_2016_2022_R200

echo %m_nam%

set input_name=static_model_inputs_no_p_o.xlsx
rem set input_name=static_model_inputs.xlsx

echo %input_name%

call activate geo_env
python 05b_plot_output_only.py %2 %m_nam% %input_name%

PAUSE