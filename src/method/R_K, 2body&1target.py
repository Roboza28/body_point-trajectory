import numpy as np
from settings import G
from src.models.solver import B2T1
from src.models.runge_kutta import rk4_step
from src.models.target_1 import f_2_center
from src.utils.utils import params_tuple


# # Н.У.

y0 =  np.array([5, 0, 0,    0., 0.322386, 0.0969581,   0, 0.48479, -1.44494])
# y0 =  np.array([5, 0, 0,    0., 1, 0.0969581,    0, 0.48479, -1.44494])
# y0 =  np.array([1/2, 0, 0,    0, 1, 0,    0, 0, 0])
# y0 = np.array([1/2, 0, -1,    0, 0.5,  0.1,    0,   0, 0]) #!! отличное решение
# y0 = np.array([1/2, 0, -1,    0, 0.1,  0.1,    0,   0, 0]) #!! отличное решение
# y0 = np.array([1/2, 0,  1,    0, 0.1,  0.5,    0,   0, 0]) #!! отличное решение
# y0 = np.array([1/2, 0, -1,    0, 0.01, 0.1,    0,   0, 0]) #!! отличное решение
# y0 = np.array([2/3, 0,  0,    0, 0.05,   1,    0.1, 0, 0]) #!! отличное решение
#y0 =  np.array([2/3, 0, 0,    0., 0.122386, 0.2269581,   0, 0.28479, -1.44494])  #!! отличное решение
#y0 =  np.array([1/2, 0, -1,    0, 0.1 , 0.1,    0, 0, 0])  #!! отличное решение
# y0 =  np.array([1/2, 0, -1,    0, 0.01 , 0.1,    0, 0, 0])  #!! отличное решение


m = 1
j = 1
b = 0.7
#b = 0.2
A = 1/2
# A = 1


yC2 = [np.array([0, 0, 0]), np.array([0.5, 0, 0])]
# yC2 = np.array([0.5, 0, 0])
#yC2 = [0,0,0]
#yC2 = [20,0,0]
#yC2 = [0.5,0,0]


mc = np.sqrt(A / G)
Rs = np.linalg.norm(yC2[1] - yC2[0])
T = 2 * np.pi * np.sqrt(Rs ** 3 / (G * mc))
print(T)

t0_task = 0
tEnd_task = 100
n_task = 1000
tau = (tEnd_task - t0_task) / n_task


params = params_tuple(m, j, b, A, yC2, 0)
solver = B2T1(f_2_center, rk4_step, t0_task, tEnd_task, tau, y0, params)
solver.solve()
solver.create_plotly_graph()
solver.create_graph_energy()
