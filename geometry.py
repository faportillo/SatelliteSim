import math
import numpy as np

'''
    Space mission geometry.
    Angles computed in radians
'''


def spherical_to_cart(r, theta, phi):
    x = r * np.sin(phi) * np.cos(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(phi)
    return x, y, z


def get_distance(vec1, vec2):
    d = np.sqrt((vec1[0] - vec2[0]) ** 2 + (vec1[1] - vec2[1]) ** 2 + (vec1[2] - vec2[2]) ** 2)
    return d


def spherical_to_vector(az, el):
    """
        NOTE: Might be deprecated
        Convert azimuth and elevation to unit vector coordinates (on the unit sphere with the craft in the center)
        :param az: azimuth - the angular distance from the north or south point of the horizon
                            to the point at which a vertical circle passing
                            through the object intersects the horizon.
        :param el: elevation - arc length distance above or below the equator
        :return: unit vector coordinates
    """
    x = math.cos(az) * math.cos(el)
    y = math.sin(az) * math.cos(el)
    z = math.sin(el)
    return x, y, z


def vector_to_spherical(x, y, z):
    """
        Convert unit vector coordinates to azimuth and elevation
        :param x: x point on unit sphere
        :param y: y point on unit sphere
        :param z: z point on unit sphere
        :return: azimuth and elevation
    """
    az = math.atan2(y, x)
    el = math.asin(z)

    return az, el


def angular_radius_spherical_earth(altitude, radius_e):
    """
        Angular radius of the spherical Earth as seen from the spacecraft (SMAD pg. 111)
        :param altitude: Altitude of spacecraft
        :param radius_e: spherical earth radius
        :return: rho - the angular radius of the spherical earth as seen from spacecraft (radians)
    """
    if radius_e is None:
        radius_e = 6378
    return math.asin(radius_e / (radius_e + altitude))


def angular_radius_center_earth(altitude, radius_e):
    """
        Angular radius measured at the center of the Earth of the region seen by the spacecraft (SMAD pg. 111)
        :param altitude: Altitude of spacecraft
        :param radius_e: spherical earth radius
        :return: lambda_0 - the angular radius measured at the center of Earth of the region seen by craft (radians)
    """
    if radius_e is None:
        radius_e = 6378
    return math.acos(radius_e / (radius_e + altitude))


def distance_horizon(radius_e, altitude=None, lambda_0=None):
    """
        Distance of the spacecraft to the Earth's horizon (SMAD pg. 111)
        :param radius_e: spherical earth radius
        :param altitude: Altitude of spacecraft
        :param lambda_0: the angular radius measured at the center of Earth of the region seen by craft
        :return: D_max - Distance of the spacecraft to the Earth's horizon
    """
    assert altitude is not None or lambda_0 is not None
    if radius_e is None:
        radius_e = 6378

    if altitude is not None:
        return ((radius_e + altitude) ** 2 - (radius_e ** 2)) ** (1 / 2)
    else:
        return radius_e * math.tan(lambda_0)


def spacecraft_to_earth_coords(altitude, radius_e, phi_e, nadir, long_ssp, lat_ssp):
    """
        Compute coordinates on the Earth of the given subsatellite point at (long_ssp, lat_ssp) (SMAD pg. 114)
         and target direction (phi_e, nadir)
        :param altitude: Altitude of spacecraft
        :param radius_e: spherical earth radius
        :param phi_e: azimuth
        :param nadir: nadir angle, measured from spacecraft ssp (subsatellite point) to target
        :param long_ssp: longitude on Earth of ssp
        :param lat_ssp: latitude on Earth of ssp
        :return: longitude and latitutde of taget on the Earth, east-west direction from ssp
    """
    if radius_e is None:
        radius_e = 6378
    rho = angular_radius_spherical_earth(altitude, radius_e)
    eps = math.acos(math.sin(nadir) / math.sin(rho))
    lambda_s = math.radians(90) - nadir - eps
    lat_p_prime = math.acos(
        math.cos(lambda_s) * math.sin(lat_ssp) + math.sin(lambda_s) * math.cos(lat_ssp) * math.cos(phi_e))
    lat_p = math.radians(90) - lat_p_prime
    delta_l = math.acos(
        (math.cos(lambda_s) - math.sin(lat_ssp) * math.sin(lat_p)) / (math.cos(lat_ssp) * math.cos(lat_p)))
    long_p = delta_l + lat_p

    return long_p, lat_p  # , 'west' if phi_e > math.radians(180) else 'east'


def earth_to_spacecraft_coords(altitude, radius_e, long_ssp, lat_ssp, long_p, lat_p):
    """

        :param altitude:
        :param radius_e:
        :param long_ssp:
        :param lat_ssp:
        :param long_p:
        :param lat_p:
        :return:
    """
    if radius_e is None:
        radius_e = 6378
    delta_l = abs(long_ssp - lat_ssp)
    rho = angular_radius_spherical_earth(altitude, radius_e)
    lambda_s = math.acos(math.sin(lat_ssp) * math.sin(lat_p) + math.cos(lat_ssp) * math.cos(lat_p) * math.cos(delta_l))
    phi_e = math.acos(
        (math.sin(lat_p) - math.cos(lambda_s) * math.sin(lat_ssp)) / (math.sin(lambda_s) * math.cos(lat_ssp)))
    n = math.atan2(math.sin(rho) * math.sin(lambda_s), (1 - math.sin(rho) * math.cos(lambda_s)))

    return phi_e, n


def subsatellite_coords_ascending(inclination, w, t, w_e=None):
    """
        Latitude and longitude relative to the ascending node for satellite in circular orbit at inclination (SMAD pg. 116)
        :param inclination: orbital inclination of satellite
        :param w: angular velocity of satellite
        :param t: time since satellite crossed the equator northbound.
        :param w_e: angular rotation of planetary body. If none, then assume it's Earth
        :return: Subsatellite Longitude adn Latitiude
    """
    if w_e is None:
        w_e = 0.000072921  # Rotation of the earth in rads/s
    lat_s = math.asin(math.sin(inclination) * math.sin(w * t))
    long_s = math.atan(math.cos(inclination) * math.tan(w * t)) - (w_e * t)

    return long_s, lat_s


def ground_track_velocity(period, radius_e=None):
    """
        Ground track velocity (SMAD pg. 116)
        :param period: Orbital period in seconds
        :param radius_e: Spherical radius of planetary body. If None, assume it's Earth
        :return: velocity
    """
    if radius_e is None:
        radius_e = 6378
    velocity_ground_track = 2 * math.pi * radius_e / period
    assert velocity_ground_track <= 7.905  # km/s

    return velocity_ground_track


def area_coverage_rate(period, l_outer=None, l_inner=None, look_one_dir=False, eps=None, rho=None):
    """
        Width of swath coverage (SMAD pg. 116-117)
        :param period: Orbital period
        :param l_outer: Effective outer horizon from ground trace
        :param l_inner: Effective inner horizon from ground trace
        :param look_one_dir: Boolean to determine if looking exclusively in
                one direction (i.e. both horizons on one size of the ground trace)
        :param eps: Elevation angle
        :param rho: Angular radius of earth
        :return: ACR - Area Coverage Rate
    """
    if eps is not None and rho is not None:
        assert l_outer is None and l_inner is None
        acr = (4 * math.pi / period) * math.cos(eps + math.asin(math.cos(eps) * math.sin(rho)))
    elif l_outer is not None and l_inner is not None:
        assert eps is None and rho is None
        if l_outer == l_inner:
            acr = (4 * math.pi / period) * math.sin(l_outer)
        else:
            acr = 2 * math.pi * ((math.sin(l_outer) + math.sin(l_inner)) if look_one_dir is False else (
                    math.sin(l_outer) + math.sin(l_inner)))
    else:
        raise Exception(
            "Invalid equation parameters. Either 'l_outer' and 'l_inner'" \
            " need valid values or 'eps' and 'rho' need valid values")

    return acr


def get_max_angles_and_max_range(rho, radius_e=None, eps_min=math.radians(5)):
    """
        Given value eps_min (typically 5 deg) and angular radius of planet w.r.t Satellite, calculate max planet central
        angle (l_max), max nadir angle, n_max, measured at sat from nadir to ground station, and max range D_max,
        which the Sat will still be in view. (SMAD pg. 119)
        :param eps_min: Minimum value of spacecraft elevation (note: like azimuth and elevation, not altitude)
        :param rho: Angular radius of planet w.r.t satellite given by radius_e / (radius_e + H)
        :param radius_e: Spherical radius of planetary body
        :return: l_max, n_max, D_max
    """
    if radius_e is None:
        radius_e = 6378
    sin_n_max = math.sin(rho) * math.cos(eps_min)  # get n_max by math.asin(sin_n)
    l_max = math.radians(90) - eps_min - math.asin(sin_n_max)
    D_max = radius_e * (math.sin(l_max) / sin_n_max)

    return l_max, math.asin(sin_n_max), D_max


def inst_orbit_pole(inclination, L_node):
    """
        Calculate instantaneous orbit pole - the pole of the orbit plane at the time of the observation given the
        inclination of the orbit and the Longitude (L_node) of the ascending node. (SMAD pg. 119)
        :param inclination: Orbital inclination
        :param L_node: Longitude of the ascending node
        :return: Longitude and Latitude of the instantaneous orbit pole
    """
    lat_pole = math.radians(90) - inclination
    long_pole = L_node - math.radians(90)

    return long_pole, lat_pole


def calc_lambda_min(long_pole, lat_pole, long_gs, lat_gs):
    """
        Calculate the minimum planetary central angle between the satellite's ground track
        and the gound station (lambda_min) (SMAD pg. 120)
        :param long_pole: Longitude of the instantaneous orbit pole
        :param lat_pole: Latitude of the instantaneous orbit pole
        :param long_gs: Longitude of the ground station
        :param lat_gs: Latitude of the ground station
        :return: sine of lambda_min
    """
    sin_l_min = math.sin(lat_pole) * math.sin(lat_gs) + math.cos(lat_pole) * math.cos(long_gs) * math.cos(
        long_gs - long_pole)

    return sin_l_min


def get_angles_and_range_closest_approach(long_pole, lat_pole, long_gs, lat_gs, rho, radius_e):
    l_min = calc_lambda_min(long_pole, lat_pole, long_gs, lat_gs)
    n_min = math.atan((math.sin(rho) * l_min) / (1 - math.sin(rho) * math.cos(math.asin(l_min))))
    eps_max = math.radians(90) - math.asin(l_min) - n_min
    d_min = radius_e * (l_min / math.sin(n_min))

    return n_min, eps_max, d_min


def max_angular_rate_sat_from_gs(d_min, v_sat, radius_e=None, altitude=None, period=None):
    """
        Get the maximum angular rate of the satellite as seen from the ground station
        :param d_min: Minimum approach distance (SMAD pg. 120)
        :param v_sat: Velocity of the satellite
        :param radius_e: Radius of planetary body
        :param altitude: Height of the satellite
        :param period: Orbital period in minutes
        :return: Max angular rate
    """
    if v_sat is not None:
        assert radius_e is None and altitude is None and period is None
        theta_max = v_sat / d_min
    elif radius_e is not None and altitude is not None and period is not None:
        assert v_sat is None
        theta_max = (2 * math.pi * (radius_e + altitude)) / (period * d_min)
    else:
        raise Exception(
            "Invalid equation parameters. Either 'v_sat' is not None or" \
            " radius_e, altitude, and period are not None.")

    return theta_max


def total_azimuth_range(lambda_min, lambda_max):
    """
        Calculate the total azimuth range that the satellite covers as seen by the ground station
        :param lambda_min: minimum planetary central angle between the satellite's ground track and the gound station
        :param lambda_max: max planetary central angle between the satellite's ground track and the gound station
        :return: total azimuth range
    """
    delta_phi = 2 * math.acos((math.tan(lambda_min)) / (math.tan(lambda_max)))

    return delta_phi


def calc_total_time_in_view(period, lambda_min, lambda_max):
    """
        Calculate the total time the satellite will be in view from the ground station (SMAD pg. 120)
        :param period: Orbital period
        :param lambda_min: minimum planetary central angle between the satellite's ground track and the ground station
        :param lambda_max: max planetary central angle between the satellite's ground track and the ground station
        :return: Total time in view (T)
    """
    T = ((period / 180) * (1 / 0.0174533)) * math.acos((math.cos(lambda_max)) / (math.cos(lambda_min)))

    return T


def calc_phi_center(lat_pole, lat_gs, lambda_min):
    """
        Calculate azimuht at the center of the viewing arc at which the elevation angle is a maximum.
        :param lat_pole: Latitude of instantaneous orbit pole
        :param lat_gs: Latitude of ground station
        :param lambda_min: minimum planetary central angle between the satellite's ground track and the gound station
        :return: Phi center
    """
    phi_pole = (math.sin(lat_pole) - math.sin(lambda_min) * math.sin(lat_gs)) / (
            math.cos(lambda_min) * math.cos(lat_gs))

    phi_center = math.radians(180) - phi_pole

    return phi_center


def max_time_in_view(period, lambda_max):
    """
        Maximum time in view that occurs when satellite passes overhead and lambda_min=0
    :param period: Orbital period in minutes
    :param lambda_max: max planetary central angle between the satellite's ground track and the ground station
    :return: T_max
    """
    T_max = period * (lambda_max / (math.radians(180)))
    return T_max
