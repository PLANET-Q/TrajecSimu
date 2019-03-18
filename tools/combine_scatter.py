import numpy as np
import os
import argparse
import glob
import openpyxl as pyxl
import traceback
import json
import simplekml
import kmlplot
import pandas as pd

parser = argparse.ArgumentParser(description='予報誤差分散による落下分散を合成して1つのkmlファイルにするスクリプト')

rootdir = os.path.dirname(os.path.abspath(__file__))

# 必須の引数
parser.add_argument('scatter_folder',
                    help='誤差分散による各地上風向風速に対する落下分散が入っているフォルダ名')
parser.add_argument('place',
                    help='射点の地名(izu_sea or noshiro_sea)')
parser.add_argument('output_dir',
                    help='出力先ディレクトリ(拡張子は省略してください)')

parser.add_argument('-f', '--forecast_point_folder',
                    help='誤差分散に用いた予報風通りに風が吹いた場合の落下分散が入っているフォルダ(プロットする場合のみ指定)')
parser.add_argument('-s', '--shift_to_point',
                    help='シフト先予報ファイル')

args = parser.parse_args()

kml_bal = simplekml.Kml()
kml_para = simplekml.Kml()
# Kmlオブジェクトに射点レギュレーションパラメータを設定
regulation_param_path = os.path.abspath(os.path.join(rootdir, '../location_parameters/', args.place+'.json'))
with open(regulation_param_path, 'r') as f:
    regulation_param = json.load(f)
rail_coord, rail_radius = kmlplot.getRailCoordByRegulation(regulation_param)
print('Launcher rail coord. of ', args.place, ': ', rail_coord)

kmlplot.setKmlByDicts(regulation_param, kml_bal)
kmlplot.setKmlByDicts(regulation_param, kml_para)

wind_direction_array = np.arange(0., 360., 22.5)
wind_speed_array = np.arange(1., 8.)
n_speeds = len(wind_speed_array)
n_directions = len(wind_direction_array)

def getScatterFromExcel(sct_file):
    sct_book = pyxl.load_workbook(sct_file)
    bal_x_sheet = sct_book['弾道 x ']
    bal_y_sheet = sct_book['弾道 y ']
    para_x_sheet = sct_book['パラ x ']
    para_y_sheet = sct_book['パラ y ']

    bal_x_array = [[float(c.value) for c in R] for R in bal_x_sheet['B2':'R8']]
    bal_y_array = [[float(c.value) for c in R] for R in bal_y_sheet['B2':'R8']]
    para_x_array = [[float(c.value) for c in R] for R in para_x_sheet['B2':'R8']]
    para_y_array = [[float(c.value) for c in R] for R in para_y_sheet['B2':'R8']]
    bal_x = np.array(bal_x_array).T
    bal_y = np.array(bal_y_array).T
    bal = np.array([bal_x, bal_y]).T
    para_x = np.array(para_x_array).T
    para_y = np.array(para_y_array).T
    para = np.array([para_x, para_y]).T
    return bal, para

def getRecordFromExcel(sct_file):
    sct_book = pyxl.load_workbook(sct_file)
    GoNogo_sheet = sct_book['Go-NoGo判定 ']
    max_alt_sheet = sct_book['最高高度 ']
    max_speed_sheet = sct_book['最高速度 ']
    max_mach_sheet = sct_book['最大マッハ数']
    launch_clear_sheet = sct_book['ランチクリア速度 ']
    para_speed_sheet = sct_book['パラ展開時速度 ']
    max_press_sheet = sct_book['最大動圧 ']

    gonogo = ['Go'] * 7
    for i, row in enumerate(GoNogo_sheet.rows):
        for cell in row:
            if cell.value == 'NoGo':
                gonogo[i] = 'NoGo'
    
    max_alt_all = [[float(c.value) for c in R] for R in max_alt_sheet['B2':'Q8']]
    max_speed_all = [[float(c.value) for c in R] for R in max_speed_sheet['B2':'Q8']]
    max_mach_all = [[float(c.value) for c in R] for R in max_mach_sheet['B2':'Q8']]
    launch_clear_all = [[float(c.value) for c in R] for R in launch_clear_sheet['B2':'Q8']]
    para_speed_all = [[float(c.value) for c in R] for R in para_speed_sheet['B2':'Q8']]
    max_press_all = [[float(c.value) for c in R] for R in max_press_sheet['B2':'Q8']]
    max_alt = np.max(max_alt_all, axis=1)
    max_speed = np.max(max_speed_all, axis=1)
    max_mach = np.max(max_mach_all, axis=1)
    launch_clear = np.min(launch_clear_all, axis=1)
    para_speed = np.max(para_speed_all, axis=1)
    max_press = np.max(max_press_all, axis=1)
    return gonogo, max_alt, max_speed, max_mach, launch_clear, para_speed, max_press


diff_bal = np.zeros((n_speeds, n_directions+1, 2))
diff_para = np.zeros((n_speeds, n_directions+1, 2))
if args.forecast_point_folder is not None:
    fore_sct_file = glob.glob(os.path.join(args.forecast_point_folder, 'output_*deg.xlsx'))
    if len(fore_sct_file) == 0:
        raise FileNotFoundError('Forecast point file was not found.')
    
    fore_sct_file = fore_sct_file[0]
    fore_bal, fore_para = getScatterFromExcel(fore_sct_file)

    if args.shift_to_point is not None:
        shift_sct_file = glob.glob(os.path.join(args.shift_to_point, 'output_*deg.xlsx'))
        if len(fore_sct_file) == 0:
            raise FileNotFoundError('Shift point file was not found.')
        
        shift_sct_file = shift_sct_file[0]
        shift_bal, shift_para = getScatterFromExcel(shift_sct_file)
        diff_bal = shift_bal - fore_bal
        diff_para = shift_para - fore_para
        for i, wind_speed in enumerate(wind_speed_array):
            color_r = int(float(i / n_speeds) * 127) + 128
            drop_coords_bal = kmlplot.dropPoint2Coordinate(shift_bal[i], rail_coord[::-1])
            drop_coords_para = kmlplot.dropPoint2Coordinate(shift_para[i], rail_coord[::-1])

            for j, fore_direction in enumerate(wind_direction_array):
                name = 'fore_' + str(wind_speed) + ' [m/s]@' + str(fore_direction) + 'deg'
                kml_bal.newpoint(name=name, coords=[drop_coords_bal[j]])
                kml_para.newpoint(name=name, coords=[drop_coords_para[j]])
    else:
        for i, wind_speed in enumerate(wind_speed_array):
            color_r = int(float(i / n_speeds) * 127) + 128
            drop_coords_bal = kmlplot.dropPoint2Coordinate(fore_bal[i], rail_coord[::-1])
            drop_coords_para = kmlplot.dropPoint2Coordinate(fore_para[i], rail_coord[::-1])

            for j, fore_direction in enumerate(wind_direction_array):
                name = 'fore_' + str(wind_speed) + ' [m/s]@' + str(fore_direction) + 'deg'
                kml_bal.newpoint(name=name, coords=[drop_coords_bal[j]])
                kml_para.newpoint(name=name, coords=[drop_coords_para[j]])


gonogo_all = []
max_alt_all = []
max_speed_all = []
max_mach_all = []
launch_clear_all = []
para_speed_all = []
max_press_all = []
#for sct_folder in scatter_folders:
for wind_direction in wind_direction_array:
    sct_file = glob.glob(os.path.join(args.scatter_folder, '{}deg/output_*deg.xlsx'.format(wind_direction)))
    if len(sct_file) == 0:
        raise FileNotFoundError('Scatter file of ' + str(wind_direction) + 'deg was not found.')
        
    sct_file = sct_file[0]
    bal, para = getScatterFromExcel(sct_file)
    bal += diff_bal
    para += diff_para

    go, alt, sp, mach, lc, v_para, prs = getRecordFromExcel(sct_file)
    gonogo_all.append(go)
    max_alt_all.append(alt)
    max_speed_all.append(sp)
    max_mach_all.append(mach)
    launch_clear_all.append(lc)
    para_speed_all.append(v_para)
    max_press_all.append(prs)

    for i, wind_speed in enumerate(wind_speed_array):
        color_r = int(float(i / n_speeds) * 127) + 128
        drop_coords_bal = kmlplot.dropPoint2Coordinate(bal[i], rail_coord[::-1])
        drop_coords_para = kmlplot.dropPoint2Coordinate(para[i], rail_coord[::-1])

        line_bal = kml_bal.newlinestring(name=(str(wind_speed)+' [m/s]@'+str(wind_direction)+'deg'))
        line_bal.style.linestyle.color = simplekml.Color.rgb(color_r, 0, 0)
        line_bal.style.linestyle.width = 2
        line_bal.coords = drop_coords_bal

        line_para = kml_para.newlinestring(name=(str(wind_speed)+' [m/s]@'+str(wind_direction)+'deg'))
        line_para.style.linestyle.color = simplekml.Color.rgb(color_r, 0, 0)
        line_para.style.linestyle.width = 2
        line_para.coords = drop_coords_para

gonogo_all = np.array(gonogo_all).T
max_alt_all = np.array(max_alt_all).T
max_speed_all = np.array(max_speed_all).T
max_mach_all = np.array(max_mach_all).T
launch_clear_all = np.array(launch_clear_all).T
para_speed_all = np.array(para_speed_all).T
max_press_all = np.array(max_press_all).T
judge_both     = pd.DataFrame(gonogo_all[:, :], index = wind_speed_array, columns = wind_direction_array[:])
max_alt        = pd.DataFrame(max_alt_all[:, :], index = wind_speed_array, columns = wind_direction_array[:])
max_vel        = pd.DataFrame(max_speed_all[:, :], index = wind_speed_array, columns = wind_direction_array[:])
max_Mach       = pd.DataFrame(max_mach_all[:, :], index = wind_speed_array, columns = wind_direction_array[:])
max_Q          = pd.DataFrame(max_press_all[:, :], index = wind_speed_array, columns = wind_direction_array[:])
v_launch_clear = pd.DataFrame(launch_clear_all[:, :], index = wind_speed_array, columns = wind_direction_array[:])
v_para_deploy  = pd.DataFrame(para_speed_all[:, :], index = wind_speed_array, columns = wind_direction_array[:])

if args.shift_to_point is None:
    excel_file = pd.ExcelWriter(os.path.join(args.output_dir, 'output.xlsx'))

# write dataframe with sheet name
judge_both.to_excel(excel_file, 'Go-NoGo判定 ')
max_alt.to_excel(excel_file, '最高高度 ')
max_vel.to_excel(excel_file, '最高速度 ')
max_Mach.to_excel(excel_file, '最大マッハ数')
v_launch_clear.to_excel(excel_file, 'ランチクリア速度 ')
v_para_deploy.to_excel(excel_file, 'パラ展開時速度 ')
max_Q.to_excel(excel_file, '最大動圧 ')
# save excel file
excel_file.save()

kml_bal.save(os.path.join(args.output_dir, 'scatter_bal.kml'))
kml_para.save(os.path.join(args.output_dir, 'scatter_para.kml'))