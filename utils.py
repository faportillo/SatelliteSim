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
    return math.asin(radius_e / (radius_e + altitude))


def angular_radius_center_earth(altitude, radius_e):
    """
        Angular radius measured at the center of the Earth of the region seen by the spacecraft (SMAD pg. 111)
        :param altitude: Altitude of spacecraft
        :param radius_e: spherical earth radius
        :return: lambda_0 - the angular radius measured at the center of Earth of the region seen by craft (radians)
    """
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
    delta_l = abs(long_ssp - lat_ssp)
    rho = angular_radius_spherical_earth(altitude, radius_e)
    lambda_s = math.acos(math.sin(lat_ssp) * math.sin(lat_p) + math.cos(lat_ssp) * math.cos(lat_p) * math.cos(delta_l))
    phi_e = math.acos(
        (math.sin(lat_p) - math.cos(lambda_s) * math.sin(lat_ssp)) / (math.sin(lambda_s) * math.cos(lat_ssp)))
    n = math.atan2(math.sin(rho) * math.sin(lambda_s), (1 - math.sin(rho) * math.cos(lambda_s)))

    return phi_e, n
