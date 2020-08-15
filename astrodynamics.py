import math
import numpy as np
from vpython import *
import geometry as geo

'''
    Kepler's Laws of Planetary Motion:
    First Law: The orbit of each planet is an ellipse, with the Sun at one focus.
    Second Law: The line joining the planet to the Sun sweeps out equal areas in equal times.
    Third Law: The square of the period of a planet is proportional to the cube of its mean distance from the Sun
'''


def law_of_gravitation(G, primary_body, orbiting_body):
    """
        Netwon's Law of Gravitation
        :param G: universal constant of gravitation
        :param primary_body: Primary orbiting body
        :param orbiting_body: Satellite or other orbiting body
        :return: return magnitude of force caused by gravity
    """
    m1 = primary_body.mass
    m2 = orbiting_body.mass
    r = orbiting_body.body.pos
    F = -G * m1 * m2 * r.hat / (mag2(r - primary_body.body.pos))

    return F


def two_body_equation(a, e, v):
    """
        Solution to the two body equation of motion for satellite orbiting planet.
        Note* may need to add other equations but right now, its the polar equation of a conic section (SMAD pg. 133)
        :param a: Semi-major axis
        :param e: Eccentricity
        :param v: Polar angle
        :return: Magnitude of position vector in terms of location in the orbit
    """
    r = a * (1 - e ** 2) / (1 + e * np.cos(v))

    return r


def energy_equation(g, m1, m2, V=None, r=None, a=None):
    """
        Calculates specific orbital energy - the sum of the kinetic energy per unit mass and potential energy
        per unit mass. Will always be a negative value. As it approaches 0, the orbit becomes parabolic.
        :param g: universal constant of gravitation
        :param m1: mass of planetary body
        :param m2: mass of satellite
        :param V: relative orbital speed
        :param r: orbital distance between bodies
        :param a: semi-major axis
        :return: Specific orbital energy
    """
    u = g * (m1 + m2)
    if V is not None and r is not None:
        assert a is None
        energy = V ** 2 / 2 - u / 2
    elif a is not None:
        assert V is None and r is None
        energy = -u / (2 * a)
    else:
        raise Exception("Invalid equation parameters. Either V and r are not None or a is not"
                        " None according to equations")

    return energy


def orbital_velocity(g, primary_body, orbiting_body, type='elliptical'):
    """
        Orbital velocity of satellite relative to planetary body.
        :param g: universal constant of gravitation
        :param primary_body: Primary orbiting body
        :param orbiting_body: Satellite or other orbiting body
        :param type: {elliptical, circular, escape}
        :return: Orbital velocity of satellite
    """
    m1 = primary_body.mass
    m2 = orbiting_body.mass
    r = np.sqrt(mag2(primary_body.body.pos - orbiting_body.body.pos))
    if orbiting_body.body_type == 'planet':
        u = g * (m1 * m2)
    else:
        u = g * m1
    if type == 'elliptical':
        assert orbiting_body.orbital_info['semimajor_axis'] is not None
        a = orbiting_body.orbital_info['semimajor_axis']
        v_o = np.sqrt(u * ((2 / r) - (1 / a)))
    elif type == 'circular':
        v_o = np.sqrt(u / r)
    elif type == 'escape':
        v_o = np.sqrt(2 * u / r)
    else:
        raise Exception("Invalid orbit type.")

    return v_o


def compute_initial_velocity(vel, inclination, arg_perigee):
    x, y, z = geo.spherical_to_cart(vel, inclination, arg_perigee)
    vel_vec = -y, x, z
    return vel_vec


def angular_momentum(pos_vector, vel_vector):
    """
        Calculate angular momentum of of satellite - the total angular momentum divided by it's mass.
    :param pos_vector: Position vector
    :param vel_vector: Velocity vector
    :return: angular momentum
    """
    h = np.cross(pos_vector, vel_vector)
    return h


def compute_semimajor_axis(u=None, eps=None, P=None, r_a=None, r_p=None):
    """
        Compute the semi-major axis, which describes the size of the ellipse.
        :param u: G * (m1 + m2)
        :param eps: specific orbital energy
        :param P: Orbital period in minutes
        :param r_a: radius of the apogee
        :param r_p: radius of the perigee
        :return: a - the semi-major axis of the orbit
    """
    if u is not None and eps is not None:
        assert r_a is None and r_p is None and P is None
        a = -u / (2 * eps)
    elif u is not None and P is not None:
        assert r_a is None and r_p is None and eps is None
        a = ((P / 2 * np.pi) ** 2 * u) ** (1 / 3)
    elif r_a is not None and r_p is not None:
        assert u is None and eps is None and P is None
        a = (r_a + r_p) / 2
    else:
        raise Exception("Invalid equation parameters.")
    return a


def compute_eccentricity_vector(v, h, r, u):
    """
        Eccentricy vector of a Kepler orbit - the dimensionless vector with direction pointing from apoapsis
        to periapsis and with the magnitude equal to the orbit's scalar eccentricity
        :param v: Velocity vector
        :param h: Specific orbital momentum vector
        :param r: position vector
        :param u: gravitation parameter - G*(m1+m2)
        :return: Eccentricity vector
    """
    e = (np.cross(v, h)) / u - (r / np.linalg.norm(r))

    return e


def compute_eccentricity(a, r_a=None, r_p=None):
    """
        Compute eccentricity of the orbit
        :param a: semi-major axis
        :param r_a: radius of the apogee
        :param r_p: radius of the perigee
        :return: eccentricity
    """
    if r_a is not None:
        assert r_p is None
        e = (r_a / a) - 1
    elif r_p is not None:
        assert r_a is None
        e = 1 - (r_p / a)
    else:
        raise Exception("Invalid equation paramets. Either r_a is not None or r_p is not None according to "
                        "equation in SMAD pg. 137")

    return e


def compute_inclination(h):
    """
        Compute inclination of the orbit (SMAD pg. 137)
        :param h: angular momentum vector
        :return: inclination
    """
    i = np.arccos((h[2] / np.linalg.norm(h)))  # Note: h[2] is the z-component of the vector

    return i


def compute_nodal_vector(h):
    """
        Compute nodal vector in the direction of the ascending node
        :param h: Angular momentum vector
        :return: nodal vector
    """
    Z = np.array([0, 0, 1])
    return np.cross(Z, h)


def compute_right_ascension(n=None, h=None):
    """
        Compute right ascension of the ascending node based on angular momentum vector(SMAD pg. 137)
        :param n:
        :return:
    """
    if n is None and h is not None:
        n = (-h[1], h[0], 0)

    omega = np.arccos((n[0] / np.linalg.norm(n))) if n[1] >= 0 else 2 * np.pi - np.arccos((n[0] / np.linalg.norm(n)))

    return omega


def compute_arg_perigee(n, e):
    """
        Argument of the perigee, the angle from the body's ascending node to its perigee in the direction of motion
        :param n: nodal vector in direction of ascending node
        :param e: eccentricity vector
        :return: Argument of perigee (w)
    """
    if e[2] != 0:
        w = np.arccos((np.dot(n, e)) / (np.linalg.norm(n) * np.linalg.norm(e)))
    else:
        w = np.arctan2(e[1], e[0])

    return w


def compute_perigee_radius(a, e):
    """
        Compute radius of perigee
        :param a: Semi-major axis
        :param e: eccentricity scalar
        :return: Radius of perigee
    """
    r_p = a * (1 - e)

    return r_p


def compute_apogee_radius(a, e):
    """
        Compute radius of apogee
        :param a: Semi-major axis
        :param e: eccentricity scalar
        :return: Radius of apogee
    """
    r_a = a * (1 + e)

    return r_a


def compute_orbital_period(a=None, u=None, long_1=None, long_2=None):
    """
        Calculate orbital period in minutes
        :param a: Semi-major axis
        :param u: Gravitation constant G*(m1+m2)
        :return: orbital period in minutes
    """
    if a is not None and u is not None:
        assert long_1 is None and long_2 is None
        P = 2 * np.pi * np.sqrt((a ** 3) / u)
    elif long_1 is not None and long_2 is not None:
        assert a is None and u is None
        delta_L = long_1 - long_2
        if delta_L > 0:  # prograde orbit
            P = 4 * (360 - delta_L)
        elif delta_L < 0:  # retrograde orbit
            P = 4 * (delta_L - 360)
        else:
            raise Exception("Not in orbit if long_1 == long_2")
    else:
        raise Exception("Invalid equation parameters")

    return P


def compute_orbital_frequency(a, u):
    """
        Compute the orbital frequency
        :param a: Semi-major axis
        :param u: Gravitation constant G*(m1+m2)
        :return: orbital frequency in rad/s
    """
    w_0 = np.sqrt(u / (a ** 3))
    return w_0


def compute_mean_motion(u, a):
    """
        Compute the mean motion - average angular velocity, determined from the semi-major axis of the orbit.
        :param u: Gravitational constant G*(m1+m2)
        :param a: semi-major axis
        :return: mean motion, n (not be confused with nodal vector)
    """
    n_ = np.sqrt(u / (a ** 3))

    return n_


def compute_true_anomaly(e, r=None, v=None, M=None):
    """
        Calculate true anomaly - the angle between the direction of the perigee and the current position of the body
        as seen from the main focus of the ellipse.
        :param r: Position vector
        :param e: Eccentricity vector
        :param v: Velocity vector
        :return: True anomaly
    """
    if r is not None and v is not None:
        assert M is None
        true_anom = np.arccos((np.dot(e, r)) / (np.linalg.norm(e) * np.linalg.norm(r)))
        if np.dot(r, v) < 0:
            true_anom = 2 * np.pi - v
    elif M is not None:
        assert r is None and v is None
        true_anom = M + 2 * e * np.sin(M) + 1.25 * (e ** 2) * np.sin(2 * M)

    else:
        raise Exception("Invalid equation parameters. See SMAD ch. 6 for more info")
    return true_anom


def compute_eccentric_anomaly(e, true_anom):
    """
        Compute the eccentric anomaly - the angular parameter that defines the position of a body that is moving
        along an elliptic Kepler orbit.
        :param e: Eccentricity vector
        :param true_anom: True Anomaly
        :return: Eccentric anomaly
    """
    E = np.arccos((e + np.cos(true_anom)) / (1 + e * np.cos(true_anom)))

    return E


def compute_mean_anomaly(e=None, E=None, M_0=None, n_=None, t=None, t_0=None):
    """
        Compute mean anomaly - the angular distance from the pericenter which a fictitious body would have if
        it moved in a circular orbit, with constant speed, in the same orbital period
        as the actual body in its elliptical orbit.
        :param e: eccentricity
        :param E: eccentric anomaly
        :param M_0: mean anomaly at time t_0
        :param n_: average angular velocity
        :param t: time
        :param t_0: initial time
        :return:
    """
    if e is not None and E is not None:
        assert M_0 is None and n_ is None and t is None and t_0 is None
        M = E - e * np.sin(E)
    elif M_0 is not None and n_ is not None and t is not None and t_0 is not None:
        assert e is None and E is None
        M = M_0 + n_ * (t - t_0)
    else:
        raise Exception("Invalid equation parameters. Refer to SMAD pg. 140 for more details")

    return M


def compute_time_of_flight(M, M_0, n_):
    """
        Compute time of flight given the current mean anomaly and the mean anomaly at t_0.
        :param M: Mean anomaly
        :param M_0: Mean anomaly at t_0
        :param n_: average angular velocity
        :return: t - t_0, the time of flight
    """
    t_t0 = (M - M_0) / n_

    return t_t0


def third_body_perturbations(i, n, tbp_type=None):
    """
        Calculate third body perturbations
        :param i: inclination of orbit in degrees
        :param n: number of orbit revolutions per day
        :param tbp_type: type of third body perturbation
        :return: Specific third body perturbation in radians/day
    """
    right_asc_tbp = {  # This is in degrees/day, need to convert to radians in final output
        'o_moon': -0.00338,
        'o_sun': -0.00154
    }
    arg_perigee_tbp = {  # This is in degrees/day, need to convert to radians in final output
        'w_moon': 0.00169,
        'w_sun': 0.00077
    }

    if tbp_type[0] == 'o':
        tbp = (right_asc_tbp[tbp_type] * np.degrees(np.cos(np.radians(i))) / n) * (1 / 57.2958)
    elif tbp_type[0] == 'w':
        tbp = (arg_perigee_tbp[tbp_type] * (4 - 5 * np.degrees(np.sin(np.radians(i)) ** 2)) / n) * (1 / 57.2958)
    else:
        raise Exception("Invalid Third Body Pertubation type")

    return tbp
