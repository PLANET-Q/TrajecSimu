# TrajecSimu

6-dof trajectory simulation for high-power rockets.  
current version: 3.0 (11/1/2018)

## 概要
Solves a 6-dof equation of motion for a trajectory of a transonic high-power rocket.  
Limited to ones without attitude/trajectory control.

Might have some problems on Windows/Linux.

## 使い方

詳細な使用方法は[こちらのページ](https://planet-q.github.io/simu-learning/)にまとめられています．

### 必要なPythonパッケージ
- numpy
- scipy
- pandas
- matplotlib
- numpy-quaternion(https://github.com/moble/quaternion)
- numba
- openpyxl
- simplekml
- pillow
- sympy

pipを用いてインストールする場合は以下のコマンドを実行して依存パッケージをインストールします

```
pip install numpy scipy pandas matplotlib numpy-quaternion numba openpyxl simplekml pillow sympy
```

### インストール

#### ZIPダウンロード
下記をクリックしてZIPダウンロード→展開

https://github.com/PLANET-Q/TrajecSimu/archive/master.zip

#### gitでインストール
```sh
$ git clone https://github.com/PLANET-Q/TrajecSimu
```

### サンプルコードを実行
cloneまたはダウンロードした`TrajecSimu`のフォルダをコマンドプロンプト/ターミナルで開き,
`python driver_sample.py`でサンプルコードを実行.

`Config_sample`内の`sample_config_camellia.csv`はロケットパラメータ設定ファイルのサンプル,
`sample_thrust_camellia.csv`はエンジン推力データ(スラストカーブ)のサンプル,
`sample_wind.csv`は数値予報風データのサンプルです。

## ライセンス

このソフトウェアはMITライセンスのもとで公開されています.

[MIT](https://github.com/PLANET-Q/TrajecSimu/blob/master/LICENSE)

## Author

[yamamotsu](https://github.com/yamamotsu)
