import numpy as np
from src.models.runge_kutta import rk4_step
from src.models.solver import B2T2, B1T2
from src.models.target_2 import f_1_center
from src.utils.utils import params_tuple


m = 1
j = 1
b = 0.7
A = 1
#B = 0.4

# yC1 = [0,0,0]
# yC2 = [0,100,0]
yC2 = np.array([0, 0, 0])
# yC2 = [np.array([0, 0, 0]), np.array([0, 100, 0])]

t0 = 0
tEnd = 3
n = 1000
tau = (tEnd-t0)/n
t = np.linspace(t0,tEnd,n)


y0 =  np.array([5, 0, 0,    0., 0.322386, 0.0969581,    0, 0.48479, -1.44494,10, 0, 0,    0., 0.322386, 0.1869581,    0, 0.48479, -1.44494])

#t, y, E = rungeKutta(f, t0, y0, tEnd, tau)
# t, y = rungeKutta(f, t0, y0, tEnd, tau)


params = params_tuple(m, j, b, A, yC2, 0)
solver = B1T2(f_1_center, rk4_step, t0, tEnd, tau, y0, params)
solver.solve()
solver.create_plotly_graph()
solver.create_graph_energy()
