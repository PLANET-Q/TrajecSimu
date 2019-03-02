import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import os
import argparse
from scipy.interpolate import interp1d
import glob
import traceback

parser = argparse.ArgumentParser(description='入力の予報風ともっともユークリッド距離の近い過去の予報風を探すスクリプト')

# 必須の引数
parser.add_argument('src_filename',
                    help='入力の予報風csvファイル名。')
parser.add_argument('search_foldername',
                    help='検索候補となる過去の予報風csvファイルが入っているフォルダ名。')


alt_axis = np.arange(300.0, 10000, 300.0)

args = parser.parse_args()

filename = str(args.src_filename).replace(' ','')
foldername = str(args.search_foldername).replace(' ', '')

search_path = os.path.join(foldername, '*')
search_filenames = glob.glob(search_path)

src_df = pd.read_csv(filename)
src_h = src_df['altitude']
src_u = src_df['Wind (from west)']
src_v = src_df['Wind (from south)']
src_w = src_df['Wind (vertical)']
src_u_func = interp1d(src_h, src_u, fill_value='extrapolate')
src_v_func = interp1d(src_h, src_v, fill_value='extrapolate')
src_w_func = interp1d(src_h, src_w, fill_value='extrapolate')
src_wind = np.c_[
            src_u_func(alt_axis),
            src_v_func(alt_axis),
            src_w_func(alt_axis)
            ]

square_sum = []
search_wind = []
for search_filename in search_filenames:
    df = pd.read_csv(search_filename)
    search_h = df['altitude']
    search_u = df['Wind (from west)']
    search_v = df['Wind (from south)']
    search_w = df['Wind (vertical)']
    u_func = interp1d(search_h, search_u, fill_value='extrapolate')
    v_func = interp1d(search_h, search_v, fill_value='extrapolate')
    w_func = interp1d(search_h, search_w, fill_value='extrapolate')
    wind = np.c_[
            u_func(alt_axis),
            v_func(alt_axis),
            w_func(alt_axis)
            ]
    search_wind.append(wind)
    
    square_sum.append(np.sum((src_wind - wind)**2))

nearest_wind_idx = np.argmin(square_sum)
print('')
print('---------------------')
print('    NEAREST WIND FILE:', search_filenames[nearest_wind_idx])
print('SUM-SQUARE LOSS VALUE:', square_sum[nearest_wind_idx])
print('---------------------')

nearest_wind = search_wind[nearest_wind_idx]

fig = plt.figure(0)
ax = fig.add_subplot(111, projection='3d')
plt.plot(nearest_wind[:,0], nearest_wind[:,1], alt_axis, label='nearest wind')
plt.plot(src_wind[:,0], src_wind[:,1], alt_axis, label='src wind')
ax.set_xlabel('U [m/s]')
ax.set_ylabel('V [m/s]')
ax.set_zlabel('altitude [m]')
plt.legend()
plt.show()