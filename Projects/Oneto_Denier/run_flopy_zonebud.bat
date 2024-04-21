rem conda init activates a conda shell which is what conda.bat does
CALL conda.bat activate geo_env
python 01e_zonebudget_diff.py
PAUSE