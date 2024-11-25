rem conda init activates a conda shell which is what conda.bat does
set VAR_1=input_write_2000_2022_R3
set VAR_2=R3_MAR_6x_diversion_for_available_flow
echo %VAR_1%
rem CALL conda.bat activate geo_env
call activate geo_env
python 03_model_connect.py %1 %VAR_1% %VAR_2%
PAUSE