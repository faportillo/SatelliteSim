import math
import numpy as np

'''
    Kepler's Laws of Planetary Motion:
    First Law: The orbit of each planet is an ellipse, with the Sun at one focus.
    Second Law: The line joining the planet to the Sun sweeps out equal areas in equal times.
    Third Law: The square of the period of a planet is proportional to the cube of its mean distance from the Sun
'''


def law_of_gravitation(G, m1, m2, r):
    """
        Netwon's Law of Gravitation
        :param G: universal constant of gravitation
        :param m1: mass of planetary body
        :param m2: mass of satellite
        :param r: distance from center of planet to satellite
        :return: return magnitude of force caused by gravity
    """
    F = -G * ((m1 * m2) / r ** 2)

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
    r = a * (1 - e ** 2) / (1 + e * math.cos(v))

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


def orbital_velocity(g, m1, m2, r, a=None, type='elliptical'):
    """
        Orbital velocity of satellite relative to planetary body.
        :param g: universal constant of gravitation
        :param m1: mass of planetary body
        :param m2: mass of satellite
        :param r: orbital distance between bodies
        :param a: semi-major axis
        :param type: {elliptical, circular, escape}
        :return: Orbital velocity of satellite
    """
    u = g * (m1 + m2)
    if type == 'elliptical':
        assert a is not None
        v = math.sqrt(u((2 / r) - (1 / a)))
    elif type == 'circular':
        v = math.sqrt(u / r)
    elif type == 'escape':
        v = math.sqrt(2 * u / r)

    return v


def angular_momentum(pos_vector, vel_vector):
    """
        Calculate angular momentum of of satellite - the total angular momentum divided by it's mass.
    :param pos_vector: Position vector
    :param vel_vector: Velocity vector
    :return: angular momentum
    """
    h = np.cross(pos_vector, vel_vector)
    return h
