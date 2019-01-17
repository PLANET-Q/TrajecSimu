import simplekml
import numpy as np

def dropPoint2Coordinate(drop_points, coord0):
    # drop_point: [x, y] distances from launch point [m]
    # coord0: [lon, lat]

    # Earth Radius [m]
    earth_radius = 6378150.0

    deg2rad = 2 * np.pi / 360.
    lat2met = deg2rad * earth_radius
    lon2met = deg2rad * earth_radius * np.cos(np.deg2rad(coord0[1]))

    drop_coords0 = np.zeros(np.shape(drop_points))
    # [lon, lat] of each drop points
    drop_coords0[:, 0] = drop_points[:, 0] / lon2met + coord0[0]
    drop_coords0[:, 1] = drop_points[:, 1] / lat2met + coord0[1]

    drop_coords = [tuple(p) for p in drop_coords0]
    return drop_coords

def getCirclePlot(center_coord, radius_meter, n_plots=100):
    theta = np.linspace(0, 2*np.pi, n_plots, endpoint=True)
    x = np.cos(theta) * radius_meter
    y = np.sin(theta) * radius_meter
    points = np.c_[x, y]
    return dropPoint2Coordinate(points, center_coord)

def output_kml(drop_points, rail_coord, wind_speeds, regulations, filename):
    # NOTE: 入力のrail_coordなどは[lat, lon]の順(TrajecSimu準拠)だが
    # kmlでは[lon, lat]の順なのでここでつかわれる関数はこの順です
    kml = simplekml.Kml()

    for regulation in regulations:
        if 'name' in regulation:
            name = regulation['name']
        else:
            name = None

        if regulation['type'] == 'circle':
            point = regulation['center'][::-1]
            kml.newpoint(
                name='center of ' + name,
                coords=[tuple(point)]
                )
            radius = regulation['radius']
            circle_pol = kml.newlinestring(name=name)
            # Linecolor: Orange
            # 0xBGR?
            circle_pol.style.linestyle.color = '0045ff'
            circle_pol.style.linestyle.width = 4

            circle_pol.coords = getCirclePlot(point, radius)

        elif regulation['type'] == 'polygon':
            points = [ tuple(p[::-1]) for p in regulation['points']]
            line = kml.newlinestring(name=name)
            line.coords = points
            # Linecolor: Yellow
            line.style.linestyle.color = '00d7ff'
            line.style.linestyle.width = 4

        elif regulation['type'] == 'line':
            point1 = tuple(regulation['point1'][::-1])
            point2 = tuple(regulation['point2'][::-1])
            line = kml.newlinestring(name=name)
            line.coords = [point1, point2]
            # Linecolor: Yellow
            line.style.linestyle.color = '00d7ff'
            line.style.linestyle.width = 4

        else:
            raise ValueError('The KML type: ' +\
                    regulation['type'] + ' is not available')

    n_speeds = len(wind_speeds)
    for i, wind_speed in enumerate(wind_speeds):
        color_r = int(float(i / n_speeds) * 127) + 128
        drop_coords = dropPoint2Coordinate(drop_points[i], rail_coord[::-1])
        line = kml.newlinestring(name=(str(wind_speed)+' [m/s]'))
        line.style.linestyle.color = simplekml.Color.rgb(color_r, 0, 0)
        line.style.linestyle.width = 2
        line.coords = drop_coords

    kml.save(filename)
