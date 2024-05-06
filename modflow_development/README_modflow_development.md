# modflow_development Directory README
The goal of the **modflow_development** directory is to hold subfolders dedicated to different MODLFOW packages and versions. The subfolders under projects will primary contain scripts used to create the input files necessary for specific packages or scripts using flopy that develop the standard regional groundwater flow models that will be used for all of the relevant Cosumnes research. and analysis while the data for these projects is maintained in the dedicated Box folder. The scripts in these subfolders may build on or utilize the scripts, modules and packages contained in the other directories such as *01_python_scripts*.

Example: 
**SFR** is a folder for all scripts relevant to the development of the stream flow routing (SFR) package using existing data sets stored in Box and creates cleaner datasets for use in flopy model development that are stored in Box as well.  
**Regional** on the other hand is a directory for creating the standard regional groundwater flow model to serve as a basis for other model development.  
