import numpy as np

from settings import G
from src.models.solver import B1T1
from src.models.runge_kutta import rk4_step
from src.models.target_1 import f_1_center
from src.utils.utils import params_tuple


# # Н.У.

# m = 1
# j = 1
# b = 0.7
# A = 1
# B = 0.4


m = q = 1
# m = 1
# q = 1

# # 1
# eta = 0.255
# b = 0.72965

# # 2
# eta = 0.1
# b = 0.1

# # 3
# j = 1
# b = 0.7

# j = b * (3 / (4 * math.sqrt(2)) + b) / m
# k2x0, k2y0, k2z0 = 0, math.sqrt(3) / 2, 0
# rx0, ry0, rz0 = 0, 0, -j * k2y0 ** 2 / (eta * q * (m * j - b ** 2))
# k1x0, k1y0, k1z0 = - k2y0 / rz0, 0, 0


# # # 3 (без учета Rmin < R < Rmax)
# # b = 0.1
b = 0.7
j = 1

# k2x0, k2y0, k2z0 = 0, 0, math.sqrt(3) / 2
# rx0, ry0, rz0 = 0, 1, 0
# k1x0, k1y0, k1z0 = math.sqrt(3) / 2, 0, 0

# k2x0, k2y0, k2z0 = 0, math.sqrt(3) / 2, math.sqrt(3) / 2
# rx0, ry0, rz0 = 1, 0, 0
# k1x0, k1y0, k1z0 = m/b*k2x0, m/b*k2y0, m/b*k2y0-0.5

k2x0, k2y0, k2z0 = 0, 0.07, np.sqrt(3) / 2
rx0, ry0, rz0 = 1, 0, 0
k1x0, k1y0, k1z0 = 0, 0.1, -0.469403


y0 = np.array([5, 0, 0, 0., 0.322386, 0.0969581, 0, 0.48479, -1.44494])
# # y0 =  np.array([5, 0, 0,    0., 0.20, 0.0969581,    0, 0.48479, -1.44494])
# # y0 =  np.array([5, 0, 0,    0., 0., 0.,    0, 0., -1.61193])
# y0 = np.array([rx0, ry0, rz0, k1x0, k1y0, k1z0, k2x0, k2y0, k2z0])

# =============================================================================

# gamma_practice = math.atan(q * (m * j - b ** 2)/(b * k2y0**2)) + 2.85*math.pi/180
# gamma_practice = math.pi/2 - math.atan(q * (m * j - b ** 2)/(b * k2y0**2))
# gamma_practice = 0
# gamma_theory = math.atan(q * (m * j - b ** 2) / (b * k2y0 ** 2))


# mc = np.sqrt(q / G)
mc = np.sqrt(1 / (q * G))
Rs = np.linalg.norm(np.array([5, 0, 0]))
T = 2 * np.pi * Rs ** (3/2) / np.sqrt(G * 2 * mc)
print(T)


t0 = 0
# tEnd = 1
tEnd = 500
n = tEnd * 100
tau = (tEnd - t0) / n


params = params_tuple(m, j, b, q, np.array([0, 0, 0]), 0)
# t, y, E, P1, P2 = runge_kutta(f, t0, y0, tEnd, tau)
# t, y, P1, P2, P3 = runge_kutta(f_1_center, t0, y0, tEnd, tau, params)
solver = B1T1(f_1_center, rk4_step, t0, tEnd, tau, y0, params)
solver.solve()
solver.create_plotly_graph()
# solver.create_graph_energy()
