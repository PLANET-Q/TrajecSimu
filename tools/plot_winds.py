import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import os
import argparse
from scipy.interpolate import interp1d

directory = '../wind_forecast_csv/izu/init21+9h'

years = np.arange(2018, 2019)
months = np.arange(3,4)
days = np.arange(9, 13)

init_time_UTC = 12
fore_time = 9

alt_axis = np.arange(0., 2400, 150.)

fig = plt.figure(0)
ax = fig.add_subplot(111, projection='3d')

wind_u_all = []
wind_v_all = []
for y in years:
    for m in months:
        for d in days:
            date_str = '{:0>4}{:0>2}{:0>2}'.format(y, m, d)
            filename = 'interp_init' + date_str + str(init_time_UTC) + '0000FH' + str(fore_time) + 'h.csv'
            filepath = os.path.join(directory, filename)
            if not os.path.exists(filepath):
                print('file ' + filepath + 'not found.')
                continue
            
            df = pd.read_csv(filepath)
            altitude = df['altitude']
            u = df['Wind (from west)']
            v = df['Wind (from south)']
            u_func = interp1d(altitude, u, fill_value='extrapolate')
            v_func = interp1d(altitude, v, fill_value='extrapolate')
            wind_u = u_func(alt_axis)
            wind_v = v_func(alt_axis)

            plt.plot(wind_u, wind_v, alt_axis, label=date_str)
            ax.set_xlabel('U [m/s]')
            ax.set_ylabel('V [m/s]')
            ax.set_zlabel('altitude [m]')

plt.legend()
plt.show()
