from decimal import Decimal
from hydrodynamics import HydroDynamics
from utils import *


class Stabilizer:
    def __init__(self, hydrodynamics: HydroDynamics):
        self.hydrodynamics = hydrodynamics
        self.time_exp_axis, self.mass_exp_axis = hydrodynamics.get_time_mass()
        assert self.mass_exp_axis is not None
        self.angles = None
        self.u_alpha = None
        self.run_params = self.hydrodynamics.run_kwargs
        self.w = None
        self.time = None
        self.time_step = None

    def approxU(self):
        pr = self.run_params.get('pr', 1)
        V0 = self.run_params.get('V0', 5)
        w = self.run_params.get('w', 10)
        dft, q_t = df(self.time_exp_axis, self.mass_exp_axis)
        # u_alpha = mass_axis[1:-1] + (V0*pr*(q_t)) / ((q_t + pr))
        u_alpha = self.mass_exp_axis[1:-1] + V0 - (pr * V0) / (q_t + pr)
        self.angles = w * dft
        self.u_alpha = u_alpha
        return self.angles, self.u_alpha

    def stableW(self, q: float, time: float, time_step: float):
        if self.angles is None:
            print("Get angles at first")
            return
        self.time_step = time_step
        dangles, du_da = df(self.angles, self.u_alpha)
        time = self.hydrodynamics.time_axis[-1] + np.arange(0, time, time_step)
        wt = np.zeros_like(time)
        intgrl = self.angles[-1]
        w = self.run_params.get('w', 10)
        for i in range(0, len(wt)):
            wt[i] = q / find_y(dangles, du_da / w, float(Decimal(intgrl) % (2 * Decimal(np.pi))))
            intgrl += wt[i] * time_step
        self.w = wt
        self.time = time
        return time, wt

    # короче нужно продолжить u_alpha до нужных чисел, сейчас оно считается по плохому промежутку
    def get_mass(self):
        new_mass = np.zeros_like(self.time)
        new_mass[0] = self.hydrodynamics.mass_axis[-1]
        time_step = self.time_step
        cur_angle = self.angles[-1]
        kwargs = self.hydrodynamics.run_kwargs
        for i in range(1, new_mass.shape[0]):
            cur_angle += self.w[i] * time_step
            u_point = find_y(self.angles, self.u_alpha, float(Decimal(cur_angle) % (2 * Decimal(np.pi))))
            print(f'u for non-stable {u_point}')
            adding_part = ((kwargs.get('V0', 5) * kwargs.get('pr', 1)) /
                                                         (new_mass[i - 1] + kwargs.get('V0', 5) -
                                                          u_point) - kwargs.get('pr', 1))
            new_mass[i] = new_mass[i - 1] + time_step * adding_part

        result_mass = np.hstack((self.hydrodynamics.mass_axis, new_mass))
        result_time = np.hstack((self.hydrodynamics.time_axis, self.time))

        return result_time, result_mass

