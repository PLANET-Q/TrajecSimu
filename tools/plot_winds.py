import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import os
import argparse
import glob
from scipy.interpolate import interp1d

parser = argparse.ArgumentParser(description='風をプロットするシンプルなスクリプト')

rootdir = os.path.dirname(os.path.abspath(__file__))

# 必須の引数
parser.add_argument('input_dir',
                    help='プロットする風プロファイルが入っているディレクトリ')
args = parser.parse_args()

alt_axis = np.arange(0., 2400, 150.)

fig = plt.figure(0)
ax = fig.add_subplot(111, projection='3d')

wind_u_all = []
wind_v_all = []

filenames = glob.glob(os.path.join(args.input_dir, '*.csv'))
if len(filenames) == 0:
    raise FileNotFoundError('No csv files were found in ' + args.input_dir)

for filename in filenames:
    df = pd.read_csv(filename)
    altitude = df['altitude']
    u = df['Wind (from west)']
    v = df['Wind (from south)']
    u_func = interp1d(altitude, u, fill_value='extrapolate')
    v_func = interp1d(altitude, v, fill_value='extrapolate')
    wind_u = u_func(alt_axis)
    wind_v = v_func(alt_axis)

    plt.plot(wind_u, wind_v, alt_axis, label=os.path.basename(filename))
    ax.set_xlabel('U [m/s]')
    ax.set_ylabel('V [m/s]')
    ax.set_zlabel('altitude [m]')

plt.legend()
plt.show()