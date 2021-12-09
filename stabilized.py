from hydrodynamics import HydroDynamics
from utils import *


class Stabilizer:
    def __init__(self, hydrodynamics: HydroDynamics):
        self.hydrodynamics = hydrodynamics
        self.time_exp_axis, self.mass_exp_axis = hydrodynamics.get_time_mass()
        assert self.mass_exp_axis is not None

    def approxU(self):
        run_params = self.hydrodynamics.run_kwargs
        pr = run_params.get('pr', 1)
        V0 = run_params.get('V0', 5)
        w = run_params.get('w', 10)
        dft, q_t = df(self.time_exp_axis, self.mass_exp_axis)
        # u_alpha = mass_axis[1:-1] + (V0*pr*(q_t)) / ((q_t + pr))
        u_alpha = self.mass_exp_axis[1:-1] + V0 - (pr * V0) / (q_t + pr)
        return w * dft, u_alpha
