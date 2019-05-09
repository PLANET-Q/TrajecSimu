import simplekml
import numpy as np
import json


def dropPoint2Coordinate(drop_points, coord0, mag_dec=-7.53):
    # drop_point: [x, y] distances from launch point [m]
    # coord0: [lon, lat]
    # mag_dec: magnetic deflection angle of lauch point(default: izu) [deg]

    # Earth Radius [m]
    earth_radius = 6378150.0

    deg2rad = 2 * np.pi / 360.
    lat2met = deg2rad * earth_radius
    lon2met = deg2rad * earth_radius * np.cos(np.deg2rad(coord0[1]))

    # Set magnetic declination
    mag_dec_rad = np.deg2rad(-mag_dec)
    mat_rot = np.array([[np.cos(mag_dec_rad), -1 * np.sin(mag_dec_rad)],
                        [np.sin(mag_dec_rad), np.cos(mag_dec_rad)]])

    drop_points_calibrated = np.dot(mat_rot, drop_points.T).T

    drop_coords0 = np.zeros(np.shape(drop_points))
    # [lon, lat] of each drop points
    drop_coords0[:, 0] = drop_points_calibrated[:, 0] / lon2met + coord0[0]
    drop_coords0[:, 1] = drop_points_calibrated[:, 1] / lat2met + coord0[1]

    drop_coords = [tuple(p) for p in drop_coords0]
    return drop_coords


def getCirclePlot(center_coord, radius_meter):
    # delta_theta = 20./radius_meter
    # n_plots = int(2*np.pi / delta_theta)
    theta = np.linspace(0, 2*np.pi, 64, endpoint=True)
    x = np.cos(theta) * radius_meter
    y = np.sin(theta) * radius_meter
    points = np.c_[x, y]
    return dropPoint2Coordinate(points, center_coord)


def setKmlCircle(
        kml,
        point_center,
        radius,
        name=None,
        plot_center=True,
        infill=False
        ):

    if plot_center:
        point = kml.newpoint(name=name, coords=[tuple(point_center)])
    else:
        point = None

    circleplots = getCirclePlot(
                    center_coord=point_center,
                    radius_meter=radius
                    )
    if infill:
        pol = kml.newpolygon(name=name)
        pol.outerboundaryis.coords = circleplots
    else:
        pol = kml.newlinestring(name=name)
        pol.coords = circleplots

    return point, pol


def setKmlByDicts(dicts, kml=None):
    if kml is None:
        kml = simplekml.Kml()

    for dict in dicts:
        if 'name' in dict:
            name = dict['name']
        else:
            name = None

        if dict['type'] == 'point':
            point = dict['coord'][::-1]
            kml.newpoint(
                name=name,
                coords=[tuple(point)]
            )
        elif dict['type'] == 'circle':
            point = dict['center'][::-1]
            radius = dict['radius']

            is_plot_center = True
            if 'hide_center' in dict:
                if dict['hide_center'] is True:
                    is_plot_center = False

            is_infill = False
            if 'infill' in dict:
                if dict['infill'] is True:
                    is_infill = True

            point, pol = setKmlCircle(
                            kml, point, radius, name=name,
                            plot_center=is_plot_center, infill=is_infill)

            pol.style.linestyle.width = 4
            pol.style.linestyle.color = '0045ff'

        elif dict['type'] == 'polygon':
            points = [tuple(p[::-1]) for p in dict['points']]
            line = kml.newlinestring(name=name)
            line.coords = points
            # Linecolor: Yellow
            line.style.linestyle.color = '00d7ff'
            line.style.linestyle.width = 4

        elif dict['type'] == 'line':
            point1 = tuple(dict['point1'][::-1])
            point2 = tuple(dict['point2'][::-1])
            line = kml.newlinestring(name=name)
            line.coords = [point1, point2]
            # Linecolor: Yellow
            line.style.linestyle.color = '00d7ff'
            line.style.linestyle.width = 4

        else:
            raise ValueError('The KML type: '+dict['type']+' is not available')

    return kml


def getRailCoordByRegulation(regulations):
    for regulation in regulations:
        if regulation['name'] == 'rail':
            return (np.array(regulation['center']), regulation['radius'])
    return None


def setKmlByJson(json_filename, kml=None, export_file=''):
    if kml is None:
        kml = simplekml.Kml()

    with open(json_filename, 'r') as f:
        dict_list = json.load(f)

    setKmlByDicts(dict_list, kml)

    if export_file is not '':
        kml.save(export_file)

    return kml


def output_kml(drop_points, rail_coord, wind_speeds, regulations, filename):
    # NOTE: 入力のrail_coordなどは[lat, lon]の順(TrajecSimu準拠)だが
    # kmlでは[lon, lat]の順なのでここでつかわれる関数はこの順です
    kml = simplekml.Kml()

    setKmlByDicts(regulations, kml)

    n_speeds = len(wind_speeds)
    for i, wind_speed in enumerate(wind_speeds):
        color_r = int(float(i / n_speeds) * 127) + 128
        drop_coords = dropPoint2Coordinate(drop_points[i], rail_coord[::-1])

        #directions = np.linspace(0, 360, len(drop_coords) - 1, endpoint=False)
        #for j, direction in enumerate(directions):
            #name = str(direction) + ' deg@' + str(wind_speed) + '[m/s]'
            #kml.newpoint(name=name, coords=[drop_coords[j]])

        line = kml.newlinestring(name=(str(wind_speed)+' [m/s]'))
        line.style.linestyle.color = simplekml.Color.rgb(color_r, 0, 0)
        line.style.linestyle.width = 2
        line.coords = drop_coords

    kml.save(filename)
