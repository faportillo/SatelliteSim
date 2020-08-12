from body import Body
import numpy as np
from satellite import Satellite
import astrodynamics as ast
import geometry as geo
from vpython import *

grid_dims = [-1000, 1000, -1000, 1000]
grid_scale = 1000000  # 1 grid unit = 1000km
G = 6.674e-11
dt = 10.1


# sun = Body("Sun", 1.989e30, 696.3, position=(-1000, 0, 0), body_type='planet', body_color=vector(1, 1, 0))

def test():
    sphere(pos=vector(0, 0, 0), radius=1, color=color.yellow)
    planet = sphere(pos=vector(10, 0, 0), radius=0.2, color=color.green, make_trail=True)
    planet.v = vector(0, 1.5, 0)
    m = 1
    G = 1
    M = 100
    dt = 0.005
    while 1:
        rate(400)
        r = planet.pos
        F = -G * M * m * r.hat / (mag2(r))
        a = F / m
        print("Force: {}, Acceleration: {}".format(F, a))

        planet.v = planet.v + a * dt
        planet.pos = planet.pos + planet.v * dt


def test2():
    G = 1
    dt = 0.0005
    # sun = Body(name="Sun", mass=100, dimensions=1, position=(0, 0, 0), body_type='planet', body_color=vector(1, 1, 0))
    earth = Body(name="Earth", mass=100, dimensions=1, position=(0, 0, 0), rotation_vel=7.292e1,
                 body_type='planet',
                 body_color=vector(1, 1, 1),
                 texture=textures.earth)
    sat_1 = Satellite(name="Sat1", mass=1, dimensions=[0.1, 0.1, 0.1],
                      orbital_info={"semimajor_axis": 1.3,
                                    "eccentricity": 0.0,
                                    "inclination": np.radians(0),
                                    "right_ascension": np.radians(0),
                                    "argument_perigee": np.radians(0)},
                      body_type="satellite",
                      body_color=color.red,
                      velocity_0=vector(0, 8.7, 0))

    sat_2 = Satellite(name="Sat2", mass=1, dimensions=[0.1, 0.1, 0.1],
                      orbital_info={"semimajor_axis": 4.3,
                                    "eccentricity": 0.0,
                                    "inclination": np.radians(20),
                                    "right_ascension": np.radians(30),
                                    "argument_perigee": np.radians(70)},
                      body_type="satellite",
                      body_color=color.green,
                      velocity_0=vector(3.427, 0.426, 1.172))

    sat_3 = Satellite(name="Sat3", mass=0.1, dimensions=[0.1, 0.1, 0.1],
                      orbital_info={"semimajor_axis": 5,
                                    "eccentricity": 0.0,
                                    "inclination": np.radians(0),
                                    "right_ascension": np.radians(0),
                                    "argument_perigee": np.radians(0)},
                      body_type="satellite",
                      body_color=color.white,
                      velocity_0=vector(0, 3, 0))

    sat_4 = Satellite(name="Sat4", mass=0.1, dimensions=[0.1, 0.1, 0.1],
                      orbital_info={"semimajor_axis": 3.2,
                                    "eccentricity": 0.0,
                                    "inclination": np.radians(-30),
                                    "right_ascension": np.radians(0),
                                    "argument_perigee": np.radians(-20)},
                      body_type="satellite",
                      body_color=color.yellow,
                      velocity_0=vector(-2.627, 4.842, 0.956))

    sat_5 = Satellite(name="Sat5", mass=0.1, dimensions=[0.1, 0.1, 0.1],
                      orbital_info={"semimajor_axis": 2.4,
                                    "eccentricity": 0.0,
                                    "inclination": np.radians(0.2),
                                    "right_ascension": np.radians(0),
                                    "argument_perigee": np.radians(-20)},
                      body_type="satellite",
                      body_color=color.magenta,
                      velocity_0=vector(6.458, 0, 0))

    sat_2_vel = ast.orbital_velocity(G, earth.mass, sat_5.mass, geo.get_distance(earth.position, sat_5.position),
                                     type='circular')
    print(sat_2_vel)

    body_list = [sat_1, sat_2, sat_3, sat_4, sat_5]
    while 1:
        rate(400)
        # r = sat_1.body.pos

        for b in body_list:
            F1 = -G * earth.mass * b.mass * b.body.pos.hat / (mag2(b.body.pos))
            a1 = F1 / b.mass
            b.body.v = b.body.v + a1 * dt
            b.body.pos = b.body.pos + b.body.v * dt

        earth.update_body(None, dt)


def universe():
    earth = Body(name="Earth", mass=5.972e24, dimensions=696340, position=(0, 0, 0), rotation_vel=7.292e-5,
                 body_type='planet',
                 body_color=vector(1, 1, 1),
                 texture=textures.earth)

    sat_1 = Satellite(name="Sat1", mass=4474, dimensions=[10000, 10000, 10000],
                      orbital_info={"semimajor_axis": 696340 + 1000,
                                    "eccentricity": 0.2,
                                    "inclination": np.radians(0),
                                    "right_ascension": np.radians(30),
                                    "argument_perigee": np.radians(20)},
                      body_type="satellite",
                      body_color=color.red,
                      velocity_0=vector(0, 7654, 0)
                      )
    print(sat_1.position)

    i = 0
    while i < 100:
        # rate(400)
        r = sat_1.body.pos
        F = -G * earth.mass * sat_1.mass * r.hat / (mag2(r))
        a = F / sat_1.mass
        sat_1.body.v = sat_1.body.v + a * dt
        sat_1.body.pos = sat_1.body.pos + sat_1.body.v * dt
        # F = ast.law_of_gravitation(G, earth.mass, sat_1.mass, ast.get_distance(earth.position, sat_1.position),
        #                          sat_1.body.pos.hat) / grid_scale
        # a = F / sat_1.mass
        print("Force: {}, Acceleration: {}".format(F, a))
        # sat_1.update_body(a, dt)
        # print(F)
        earth.update_body(None, dt)
        i += 1


def main():
    # test()
    test2()
    # universe()


if __name__ == '__main__':
    main()
