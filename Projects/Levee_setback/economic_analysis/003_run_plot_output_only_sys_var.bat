rem conda init activates a conda shell which is what conda.bat does


rem version with no p_o and 2019 variables
set m_nam=input_write_2014_2025_R300
rem set m_nam=input_write_2014_2025_R304
set m_nam=input_write_2014_2025_R303


echo %m_nam%

set input_name=static_model_inputs_no_p_o.xlsx

echo %input_name%

call activate geo_env
python 05b_plot_output_only.py %2 %m_nam% %input_name%

PAUSE