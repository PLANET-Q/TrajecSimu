import numpy as np
import numpy.linalg as LA
import sympy
# direction_std = 0.5 * np.pi


def getConditionalBivariateNorm(mu_X, mu_Y, sigma_XX, sigma_XY, sigma_YY, X):
    # Y|X が多変量正規分布N(mu, Sigma)に従う場合のmu, Sigmaを求める
    sigma_XX_i = LA.inv(sigma_XX)
    mu = mu_Y + np.dot(np.dot(sigma_XY.T, sigma_XX_i), (X - mu_X))
    Sigma = sigma_YY - np.dot(np.dot(sigma_XY.T, sigma_XX_i), sigma_XY)
    return mu, Sigma


def getEllipseParameters(mu, sigma, alpha=0.95):
    # mu: 平均 ([alt][x|y])
    # sigma: 分散共分散行列
    # alpha: 楕円内に点が入る確率

    sigma_i = LA.inv(sigma)
    eigvals, M = LA.eig(sigma_i)
    lmda2 = -2*np.log(1.00000000 - alpha)
    scale = np.sqrt(lmda2/eigvals)
    return scale, M, mu


def getEllipsePlot(a, b, Mrot=None, angle_rad=0, x0=0, y0=0, n_plots=100):
    theta = np.linspace(0, 2 * np.pi, n_plots)

    v0 = np.array([x0, y0])
    v = np.array([a * np.cos(theta), b * np.sin(theta)])

    if Mrot is None:
        Mrot = np.array([
            [np.cos(angle_rad), -np.sin(angle_rad)],
            [np.sin(angle_rad), np.cos(angle_rad)]
            ])

    return np.dot(Mrot, v) + v0[:, None]


def getAzimuthWindByPlot(U, V, azimuth_rad):
    # 基準高度の楕円上の各風の方位計算[-pi, +pi]
    wind_dir = np.arctan2(-U, -V)
    idx = np.argmin(abs(wind_dir - azimuth_rad))
    return np.array([U[idx], V[idx]])


def getAzimuthWind(center_pos, scale_vec, rotateMatrix, wind_direction_rad):
    unit_vec = np.array([np.sin(wind_direction_rad),
                         np.cos(wind_direction_rad)])
    rotated_unit_vec = np.dot(rotateMatrix.T, unit_vec)
    return np.dot(rotateMatrix, rotated_unit_vec * scale_vec) + center_pos


def getProbEllipse(mu, sigma, alpha=0.95, n_plots=100):
    # mu: 平均 ([alt][x|y])
    # sigma: 分散共分散行列
    # alpha: 楕円内に点が入る確率
    # n_plots: プロット数

    scale, M, _ = getEllipseParameters(mu, sigma, alpha)
    plots = getEllipsePlot(
        scale[0],
        scale[1],
        M,
        x0=mu[0],
        y0=mu[1],
        n_plots=n_plots)
    return plots[0], plots[1]


def getStatWindVector(
        statistics_parameters,
        wind_direction_deg,
        wind_std=None
        ):
    alt_axis = statistics_parameters['alt_axis']
    # alt_std = alt_axis[alt_index_std]
    n_alt = len(alt_axis)

    mu4 = np.array(statistics_parameters['mu4'])
    sigma4 = np.array(statistics_parameters['sigma4'])

    if wind_std is None:
        # ----------------------------
        # 基準風の導出
        # ----------------------------
        alt_index_std = statistics_parameters['altitude_idx_std']
        mu_stdalt = mu4[alt_index_std][2:]
        sigma_stdalt = sigma4[alt_index_std][2:, 2:]
        scale, M, _ = getEllipseParameters(mu_stdalt, sigma_stdalt, alpha=0.95)
        # wind_std_tmp = getAzimuthWindByEllipse(
        #                     mu_stdalt,
        #                     scale,
        #                     M,
        #                     np.deg2rad(wind_direction_deg)
        #                     )
        wind_std_tmp = getAzimuthWind(mu_stdalt, scale, M, np.deg2rad(wind_direction_deg))
        wind_std = np.broadcast_to(wind_std_tmp, (n_alt, 2))

    # ----------------------------
    # Probabillity Ellipse
    # ----------------------------
    stat_wind_u = []
    stat_wind_v = []
    # print('X: ', wind_std)
    for h in range(n_alt):

        # u,vが決まった時のdu,dvの条件付き正規分布の平均と共分散行列
        mu, sigma = getConditionalBivariateNorm(
            mu_X=mu4[h][2:],
            mu_Y=mu4[h][:2],
            sigma_XX=sigma4[h][2:, 2:],
            sigma_XY=sigma4[h][2:, :2],
            sigma_YY=sigma4[h][:2, :2],
            X=wind_std[h])

        scale, M, _ = getEllipseParameters(mu, sigma, alpha=0.95)
        # print('wind direction deg', wind_direction_deg)
        # w = getAzimuthWindByEllipse(
        #                     mu + wind_std[h],
        #                     scale,
        #                     M,
        #                     np.deg2rad(wind_direction_deg),
        #                     mu + wind_std[h]
        #                     )
        w = getAzimuthWind(
            mu + wind_std[h],
            scale,
            M,
            np.deg2rad(wind_direction_deg))

        stat_wind_u.append(w[0])
        stat_wind_v.append(w[1])
        print(
            'Altitude: ', alt_axis[h], ' ',
            stat_wind_u[h], ', ', stat_wind_v[h])

    stat_wind = np.array([stat_wind_u, stat_wind_v])
    return stat_wind
