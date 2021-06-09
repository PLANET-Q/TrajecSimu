import matplotlib.pyplot as plt
import numpy as np
import argparse
from scipy import fftpack
import glob
import traceback

parser = argparse.ArgumentParser(description='スラストカーブを表示するスクリプトです。')

# 必須の引数
parser.add_argument('filename',
                    help='スラストカーブのcsvファイル名。","で区切るか"*"を指定すると該当する全てのスラストカーブを重ねて表示します')
# 任意引数
parser.add_argument('-t', '--sampling_dt', help='サンプリング間隔[s] 指定しない場合スラストカーブから推測されます。')
parser.add_argument('-f', '--lpf_freq', help='LPFの遮断周波数. 指定しない場合はLPFを使用しない。')
parser.add_argument('-o', '--fitting_order', help='多項式近似の次数. 指定しない場合は多項式近似を使用しない。')

args = parser.parse_args()

filenames = str(args.filename).replace(' ','').split(',')

filename_array = []
time_array = []
thrust_raw = []
thrust_lpf = []
thrust_fit = []
for filename in filenames:
    # filenameに'*'が含まれる場合マッチするファイルが複数になる
    match_files = glob.glob(filename)

    print('--------------------------------')
    print(len(match_files), ' files matched by ' + filename)
    print('--------------------------------')
    for matchfile in match_files:
        print('  ', matchfile, ': ', end='')

        try:
            rawdata = np.loadtxt(matchfile, comments='#', delimiter=',')
        except:
            traceback.print_exc()
            continue

        filename_array.append(matchfile)

        rawtime = rawdata[:, 0]
        rawthrust = rawdata[:, 1]
        n_samples = len(rawtime)

        print('n_samples=', n_samples, ', ', end='')
        thrust_raw.append(rawthrust)
        time_array.append(rawtime)

        if args.sampling_dt is None:
            thrust_dt = rawtime[1] - rawtime[0]
        else:
            thrust_dt = float(args.t)

        print('dt=', thrust_dt, '[s], ', end='')

        if args.lpf_freq is not None:
            f_cutoff = int(args.lpf_freq)
            tf = fftpack.fft(rawthrust)
            freq = fftpack.fftfreq(n=n_samples, d=thrust_dt)
            tf2 = np.copy(tf)
            tf2[np.abs(freq) > f_cutoff] = 0
            lpf = np.real(fftpack.ifft(tf2))
            thrust_lpf.append(lpf)

        if args.fitting_order is not None:
            fit = np.polyfit(rawtime, rawthrust, int(args.fitting_order))
            func = np.poly1d(fit)

            thrust_fit.append(func(rawtime))

        print('')


for i in range(len(thrust_raw)):
    plt.figure(i)
    plt.plot(time_array[i], thrust_raw[i], label='raw')
    if args.lpf_freq is not None:
        plt.plot(time_array[i], thrust_lpf[i], label='LPF')

    if args.fitting_order is not None:
        plt.plot(time_array[i], thrust_fit[i], label='Fitted')

    plt.grid()
    plt.legend()
    plt.title(filename_array[i])

plt.show()
