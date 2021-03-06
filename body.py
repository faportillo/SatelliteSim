import numpy as np
import astrodynamics as ast
from vpython import *


class Body(object):
    def __init__(self, name, mass, dimensions, position, rotation_vel=None, orbital_info=None, body_type='planet',
                 velocity_0=None, body_color=None, texture=None):
        """
            Initialize the Body class. Will be assigned initial position at the ascending node based on orbital info.
            :param name: string name of Body
            :param mass: the body's mass in kg
            :param dimensions: the body's dimensions as a list of max size 3. If size is 1, spherical object. Size of 2
                will be oblong sphere object, size of 3 will be cube object (typically for satellites)
            :param position: The initial position
            :param orbital_info: Dict with orbital characteristics {semimajor_axis, eccentricity, inclination,
                right_ascension, argument_perigee}. *Note: rotation_vel in radians/sec
            :param body_type: Type of body {'planet', 'satellite'}
            :param velocity_0: Initial velocity vector of the body, defaults to None
            :param altitude_0: The initial altitude of the body
            :param body_color: The body's color, random if not specified
            :param texture: Can set body to have a texture instead of a color
        """
        self.name = name
        self.mass = mass
        self.dimensions = dimensions
        self.rotation_vel = rotation_vel
        self.orbital_info = orbital_info
        self.body_type = body_type
        self.body_color = body_color
        self.texture = texture

        if self.orbital_info is not None:
            print(self.orbital_info)
            for key in self.orbital_info:
                print(self.name, "->", key, ":", self.orbital_info[key])

        self.position = position

        print("{} -> Initial Position: {}".format(self.name, self.position))
        if body_type == 'planet':
            self.body = sphere(pos=vector(self.position[0], self.position[1], self.position[2]),
                               radius=self.dimensions,
                               color=body_color if body_color is not None else vector.random(),
                               make_trail=True,
                               texture=texture)
        self.body.v = velocity_0

    def update_body(self, a, dt):
        """
            Update the position and velocity of the body.
            :param a: Acceleration, based on Newton's Law of Gravitation
            :param velocity_vector: velocity vector calculated from other numerical methods
            :param dt: Global time interval
            :return: None, just updates position and velocity in VPython
        """
        if a is not None:
            self.body.v = self.body.v + a * dt
            self.body.pos = self.body.pos + self.body.v * dt

        if self.body_type == 'planet':
            self.body.rotate(angle=radians(self.rotation_vel * dt), axis=vector(0, 1, 0))
