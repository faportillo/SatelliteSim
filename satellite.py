import numpy as np
import math
import astrodynamics as ast
import geometry as geo
from vpython import *
from body import Body


class Satellite(Body):
    def __init__(self, name, mass, dimensions, position=None, rotation_vel=None, orbital_info=None,
                 body_type='satellite', orbit_type='circular', primary_body=None, G=6.674e-11, body_color=None,
                 texture=None):
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
        self.name = name
        self.mass = mass
        self.dimensions = dimensions
        self.rotation_vel = rotation_vel
        self.body_type = body_type
        self.body_color = body_color
        self.texture = texture

        if position is not None:
            self.position = position
        elif position is None and orbital_info is not None:
            self.orbital_info = orbital_info
            self.r_p = ast.compute_perigee_radius(self.orbital_info["semimajor_axis"],
                                                  self.orbital_info["eccentricity"])
            print("{} -> Radius of Perigee: {}".format(name, self.r_p))
            # self.r_p = ast.two_body_equation(self.orbital_info["semimajor_axis"], self.orbital_info["eccentricity"],
            #                                 self.orbital_info["argument_perigee"] + self.orbital_info[
            #                                     "right_ascension"])
            self.position = geo.cart_to_vpython(geo.spherical_to_cart(self.r_p, self.orbital_info["inclination"],
                                                                      self.orbital_info['argument_perigee']))

        self.body = box(pos=vector(self.position[0], self.position[1], self.position[2]),
                        length=dimensions[0], height=dimensions[1], width=dimensions[2],
                        color=body_color if body_color is not None else vector.random(),
                        make_trail=True,
                        texture=texture)

        mag_velocity = ast.orbital_velocity(G, primary_body, self, type=orbit_type)
        # velocity_vec = geo.cart_to_vpython(geo.spherical_to_cart(mag_velocity, self.orbital_info['inclination'],
        #                                                         self.orbital_info['argument_perigee']))
        velocity_vec = self.body.pos.hat * mag_velocity
        self.velocity_0 = velocity_vec  # vector(velocity_vec[0], velocity_vec[1], velocity_vec[2])
        print("{} -> Initial Velocity: {}".format(self.name, self.velocity_0))

        super().__init__(name=self.name, mass=self.mass, dimensions=self.dimensions, position=self.position,
                         rotation_vel=self.rotation_vel, orbital_info=self.orbital_info, body_type=self.body_type,
                         velocity_0=self.velocity_0, body_color=self.body_color, texture=self.texture)
