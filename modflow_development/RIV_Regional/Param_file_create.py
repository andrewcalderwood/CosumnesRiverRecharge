# -*- coding: utf-8 -*-
"""
Created on Sat Nov 28 19:17:36 2020

@author: arodri44
"""

import numpy as np

# Read in the file
File = open('Real1.PAR', 'r')
list_of_lines = File.readlines()
list_of_lines [2] = "Real2.DAT\n"

File = open('Real1.PAR','w')
File.writelines(list_of_lines)
File.close()