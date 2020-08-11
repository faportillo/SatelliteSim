import numpy as np
from vpython import *


class Body(object):
    def __init__(self, name, mass, dimensions, position=None, orbital_info=None, body_type='planet', velocity_0=None,
                 altitude_0=1000, body_color=None, texture=None):
        """
            Initialize the Body class. Will be assigned initial position at the ascending node based on orbital info.
            :param name: string name of Body
            :param mass: the body's mass in kg
            :param dimensions: the body's dimensions as a list of max size 3. If size is 1, spherical object. Size of 2
                will be oblong sphere object, size of 3 will be cube object (typically for satellites)
            :param position: The initial position
            :param orbital_info: Dict with orbital characteristics {semimajor_axis, eccentricity, inclination,
                right_ascension, argument_perigee}
            :param body_type: Type of body {'planet', 'satellite'}
            :param velocity_0: Initial velocity vector of the body, defaults to None
            :param altitude_0: The initial altitude of the body
            :param body_color: The body's color, random if not specified
            :param texture: Can set body to have a texture instead of a color
        """
        self.name = name
        self.mass = mass
        self.dimensions = dimensions
        self.orbital_info = orbital_info
        self.altitude_0 = altitude_0
        if position is not None:
            self.position = position
        elif orbital_info is not None:  # If initial position not specified, then assume it's orbiting body
            self.position = (self.altitude_0 * np.cos(self.orbital_info['right_ascension']),
                             self.altitude_0 * np.cos(self.orbital_info['right_ascension']), 0)

        self.velocity = velocity_0

        if body_type == 'planet':
            self.body = sphere(pos=self.position,
                               size=self.dimensions,
                               color=body_color if body_color is not None else vector.random(),
                               make_trail=True,
                               texture=texture)
        elif body_type == 'satellite':
            self.body = box(pos=self.position,
                            length=self.dimensions[0], height=self.dimensions[1], width=self.dimensions[2],
                            color=body_color if body_color is not None else vector.random(),
                            make_trail=True,
                            texture=texture)

    def update_body(self, a, velocity_vector, dt):
        """
            Update the position and velocity of the body.
            :param a: Acceleration, based on Newton's Law of Gravitation
            :param velocity_vector: velocity vector calculated from other numerical methods
            :param dt: Global time interval
            :return: None, just updates position and velocity in VPython
        """
        self.velocity = velocity_vector + a * dt
        self.position = self.position + velocity_vector * dt
