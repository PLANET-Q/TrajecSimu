#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 17:33:18 2017

@author: shugo
"""

# sample driver script for trajectory simulation

import numpy as np
from Rocket_simu import Rocket_simu

# define path and filename of a csv file
# csv_filename = 'Parameters_csv/2018izu/relenza_ntk.csv'
# csv_filename = 'Parameters_csv/newmitei_parameters_C63.csv'
# csv_filename = 'Parameters_csv/felix_2018izu_parameters.csv'
csv_filename = 'Parameters_csv/2018noshiro/felix_2018noshiro.csv'

# create an instance
myrocket = Rocket_simu(csv_filename)

# ------------------------------------
# run a single trajectory computation 
# ------------------------------------
myrocket.run_single()

# ------------------------------------
# run a loop for landing point distribution
# ------------------------------------
#"""
wind_direction_array = np.linspace(0.,360.,9)  # wind direction array
wind_speed_array = np.linspace(1.,7., 7)       # wind speed array
myrocket.run_loop(wind_direction_array, wind_speed_array)
#"""






