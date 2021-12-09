import numpy as np


def generate_array(a1: float, a2: float):
    angles_36 = np.array([float(i) for i in range(0, 36)])
    y = np.array([])
    for i in range(0, 50):
        y = np.concatenate((y, (y[-1] if y.size != 0 else 0) + a1 * angles_36))
        y = np.concatenate((y, y[-1] + a2 * angles_36))
    return y


def generate_smooth(a1: float, a2: float, angles):
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


def df(x: np.array, y: np.array):
    assert x.size == y.size
    return np.array(x[1:-1]), np.array([(y2-y0)/(x2-x0) for x2, x0, y2, y0 in zip(x[2:], x, y[2:], y)])


