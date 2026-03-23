import numpy as np
from src.models.runge_kutta import rk4_step
from src.models.solver import B2T1_matpoint
from src.models.target_1 import f_2_center_mat_point
from src.utils.utils import params_tuple


m = 1
A = 1/2

x0 = 1
y0 = 0
z0 = 0
vx0 = 0
vy0 = 1.1
vz0 = 0

y0 = np.array([x0, vx0, y0, vy0, z0, vz0])
yC2 = [np.array([0, 0, 0]), np.array([0.1, 0, 0])]

t0 = 0
tEnd = 100
n = tEnd * 100
tau = (tEnd - t0) / n

params = params_tuple(m, 0, 0, A, yC2, 0)
solver = B2T1_matpoint(f_2_center_mat_point, rk4_step, t0, tEnd, tau, y0, params)
solver.solve()
solver.create_plotly_graph()
solver.create_graph_energy()
