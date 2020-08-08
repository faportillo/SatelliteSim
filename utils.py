import numpy as np
import math


def spherical_to_vector(az, el):
    """
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
        :return: Subsatellite Latitiude and Longitude
    """
    if w_e is None:
        w_e = 0.000072921  # Rotation of the earth in rads/s
    lat_s = math.asin(math.sin(inclination) * math.sin(w * t))
    long_s = math.atan(math.cos(inclination) * math.tan(w * t)) - (w_e * t)

    return lat_s, long_s


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


def get_max_angles_and_range(rho, radius_e=None, eps_min=math.radians(5)):
    """
        Given value eps_min (typically 5 deg) and angular radius of planet w.r.t Satellite, calculate max planet central
        angle (l_max), max nadir angle, n_max, measured at sat from nadir to ground station, and max range D_max,
        which the Sat will still be in view.
    :param eps_min: Minimum value of spacecraft elevation (note: like azimuth and elevation, not altitude)
    :param rho: Angular radius of planet w.r.t satellite given by radius_e / (radius_e + H)
    :param radius_e: Spherical radius of planetary body
    :return: l_max, n_max, D_max
    """
    if radius_e is None:
        radius_e = 6378
    sin_n = math.sin(rho) * math.cos(eps_min)  # get n_max by math.asin(sin_n)
    l_max = math.radians(90) - eps_min - math.asin(sin_n)
    D_max = radius_e * (math.sin(l_max) / sin_n)
    
    return l_max, math.asin(sin_n), D_max
