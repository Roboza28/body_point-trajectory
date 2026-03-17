import numpy as np

from src.models.runge_kutta import rk4_step
from src.models.solver import B1T1
from src.models.target_1 import f_1_center_and_magnetic_field
from src.utils.utils import params_tuple


# # Н.У.

y0 =  np.array([5, 0, 0,    0., 0.322386, 0.0969581,    0, 0.48479, -1.44494])


m = 1
j = 1
b = 0.7
A = 1
B = 0.1

t0 = 0
tEnd = 500
n = tEnd* 100
tau = (tEnd-t0)/n
t = np.linspace(t0,tEnd,n)


params = params_tuple(m, j, b, A, [], B)
solver = B1T1(f_1_center_and_magnetic_field, rk4_step, t0, tEnd, tau, y0, params)
solver.solve()
solver.create_plotly_graph()
