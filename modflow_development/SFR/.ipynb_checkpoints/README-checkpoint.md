# SFR script folder README

- SFR_reach_input_setup is the original script for SFR reach input that samples 10 m DEM elevation directly at NHD line locations sampled every 10 m.
- SFR_XS_to_reach is the newest version of script to identify the minimum reach elevations in a more fluid way. This works for both the Cosumnes River and Deer Creek. This script supercedes SFR_reach_input_setup.py which used a simpler elevation sampling directly along the NHD lines.
- SFR_XS_setup was the original testing ground for auto cross-sections, and helps identify the spatial location of Constantine's cross-sections, also creates Michigan Bar rating curve table
- sfr_grid_comparison plots the profile of the sfr grid cells between the local Oneto-Denier and regional model which use the 2m and 10 m DEMs respectively. The local model shows more channel incision (1-2m) than the regional model predicts.
