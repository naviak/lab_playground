import numpy as np
import matplotlib.pyplot as plt

time_step = 0.01
time = 30
a1 = 1.1
a2 = 0.7
V0 = 5
w = 10
pr = 2
init_mass = 0.

# a1 * a + a2 * sina
# a1 > a2 строго моннотонна
# в функцию потока

def generate_array():
    angles_36 = np.array([float(i) for i in range(0,36)])
    y = np.array([])
    for i in range(0,50):
        y = np.concatenate((y , (y[-1] if y.size != 0 else 0) + a1 * angles_36))
        y = np.concatenate((y , y[-1] + a2 * angles_36)) 
    return y


def generate_smooth(angles):
    y = np.zeros(angles.shape[0])
    for i in range(1, len(y)):
        y[i] = y[-1] + a1 * angles[i] + a2 * np.sin(angles[i])
    return y


def find_y(x, y, x_point):
    dists = np.abs(x - x_point)
    inds = dists.argsort()[:2]
    closest_ys = np.sort(y[inds])
    closest_xs = np.sort(x[inds])
    val = (x_point - closest_xs[0]) / (closest_xs[1] - closest_xs[0]) * (closest_ys[1] - closest_ys[0]) + closest_ys[0]
    return val


def get_mass(time_axis, angles, q):
    global time_step, init_mass, pr, V0, w
    mass_axis = np.zeros_like(time_axis)
    mass_axis[0] = init_mass
    for i in range(1, mass_axis.shape[0]):
        angle = i * time_step * w 
        u = find_y(angles, q, angle)
        print(f'u for nonstable {u}')
        mass_axis[i] = mass_axis[i-1] + time_step * ((V0 * pr) / (mass_axis[i-1] + V0 - u) - pr)
    return mass_axis


def df(x: np.array, y: np.array):
    assert x.size == y.size
    return np.array(x[1:-1]), np.array([(y2-y0)/(x2-x0) for x2, x0, y2, y0 in zip(x[2:], x, y[2:], y)])

def get_new_u(time_axis, mass_axis):
    dft, q_t = df(time_axis,mass_axis) 
    u_alpha = mass_axis[1:-1] + (V0*pr*(q_t)) / ((q_t + pr))
    return w * dft, u_alpha

def get_new_w(angles, u_alpha, q):
    global time_step, w
    dangles, du_da = df(angles, u_alpha)
    #plt.plot(dangles, du_da)
    wt = np.zeros_like(dangles)
    dt = time_step
    wt[0] = 10
    intgrl = wt[0] * dt
    for i in range(1, len(wt)):
        intgrl += wt[i] * dt
        wt[i] = 25 / find_y(dangles, du_da, intgrl)
    return dangles , wt

def get_mass_mod(time_axis, angles, q, wt):
    global time_step, init_mass, pr, V0, w
    mass_axis = np.zeros_like(time_axis)
    mass_axis[0] = init_mass
    cur_angle = 0
    for i in range(1, mass_axis.shape[0]):
        angle = time_step * wt[i]
        cur_angle += angle
        u = find_y(angles, q, cur_angle)
        print(f'u for stable {u}')
        adding_part = time_step * ((V0 * pr) / (mass_axis[i-1] + V0 - u) - pr)
        mass_axis[i] = mass_axis[i-1] + adding_part
    return mass_axis

angles = np.arange(0,720, 0.01)
u_exp = generate_smooth(angles)
dangle, du_a = df(angles, u_exp)
#plt.plot(dangle, du_a) # экспериментальная производная
time_axis = np.arange(0, time, time_step)

mass_axis = get_mass(time_axis, angles, u_exp)

new_angles, u_alpha = get_new_u(time_axis, mass_axis)

angs, wt = get_new_w(new_angles, u_alpha,0.1)
#plt.plot(angs, wt)

time_axis_mod = np.arange(0, time // 2, time_step)
mass_axis_mod = get_mass_mod(time_axis_mod, angles, u_exp, wt)
#plt.plot(time_axis_mod, mass_axis_mod)
#plt.plot(time_axis, mass_axis)
plt.plot(angles, u_exp)
plt.plot(new_angles, u_alpha)

plt.show()