import numpy as np
import matplotlib.pyplot as plt
from utils import *
from hydrodynamics import HydroDynamics
from stabilized import Stabilizer

if __name__ == '__main__':
    w0 = 10
    t = 2 * np.pi / w0

    pump = HydroDynamics(4 * t, 0.001, 8 * np.pi, 0.001)
    pump.setGeometricU(a1=1.1, a2=0.7)

    run_params = {'V0': 5, 'init_mass': 1., 'pr': 1, 'w': w0}
    pump.run(**run_params)

    time_axis, mass_axis = pump.get_time_mass()
    exp_angles, u_exp = pump.get_angles_u()

    stabilizedPump = Stabilizer(pump)
    angles, u_alpha = stabilizedPump.approxU()
    fig, ax = plt.subplots()
    plt.plot(angles, u_alpha, label="Экспериментальная функция геометрии")
    plt.plot(exp_angles, u_exp, label="Восстановленная функция геометрии")
    ax.legend()
    fig.set_figheight(5)
    fig.set_figwidth(8)

    plt.show()
