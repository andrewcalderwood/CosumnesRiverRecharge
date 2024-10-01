rem conda init activates a conda shell which is what conda.bat does
CALL conda.bat activate geo_env
python multiprocess_recharge.py
PAUSE