call activate geo_env
jupytext --to notebook "01b_Flow and head obs plotting.py" --execute
rem jupyter nbconvert --to notebook --execute check.ipynb
jupyter nbconvert --to html "01b_Flow and head obs plotting.ipynb" --no-input

python 30_update_output.py
PAUSE