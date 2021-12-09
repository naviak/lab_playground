import numpy as np
from utils import *


class HydroDynamics:
    def __init__(self, time_of_experiment: float, time_step: float, angles: float, angles_step: float):

        self.time_axis = np.arange(0, time_of_experiment, time_step)
        self.angles = np.arange(0, angles, angles_step)
        self.mass_axis = np.zeros_like(self.time_axis)
        self.u = None
        self.run_kwargs = None

    def setGeometricU(self, a1, a2, generator=generate_smooth):
        self.u = generator(a1, a2, self.angles)

    def run(self, **kwargs):
        self.run_kwargs = kwargs
        self.mass_axis[0] = kwargs.get('init_mass', 0)
        time_step = self.time_axis[1] - self.time_axis[0]
        for i in range(1, self.mass_axis.shape[0]):
            angle = i * time_step * kwargs.get('w', 10)
            u_point = find_y(self.angles, self.u, angle)
            print(f'u for non-stable {u_point}')
            self.mass_axis[i] = self.mass_axis[i - 1] + time_step * ((kwargs.get('V0', 5) * kwargs.get('pr', 1)) /
                                                                     (self.mass_axis[i - 1] + kwargs.get('V0', 5) -
                                                                      u_point) - kwargs.get('pr', 1))

    def get_time_mass(self):
        return self.time_axis, self.mass_axis

    def get_angles_u(self):
        return self.angles, self.u
