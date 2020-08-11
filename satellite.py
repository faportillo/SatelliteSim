import numpy as np
import math
from body import Body


class Satellite(Body):
    def __init__(self, altitude):
        self.altitude = altitude
