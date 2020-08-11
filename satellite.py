import numpy as np
import math
from body import Body


class Satellite(Body):
    def __init__(self, name, mass, dimensions, position=None, orbital_info=None, body_type='planet', velocity_0=None,
                 altitude_0=1000, body_color=None, texture=None):
        """
            Satellite class
            :param name: string name of Body
            :param mass: the body's mass in kg
            :param dimensions: the body's dimensions as a list of max size 3. If size is 1, spherical object. Size of 2
                will be oblong sphere object, size of 3 will be cube object (typically for satellites)
            :param orbital_info: Dict with orbital characteristics {semimajor_axis, eccentricity, inclination,
                right_ascension, argument_perigee}
            :param body_type: Type of body {'planet', 'satellite'}
            :param velocity_0: Initial velocity vector of the body, defaults to None
            :param altitude_0: The initial altitude of the body
            :param body_color: The body's color, random if not specified
            :param texture: Can set body to have a texture instead of a color
        """
        super().__init__(name, mass, dimensions, position, orbital_info, body_type, velocity_0,
                         altitude_0, body_color, texture)
