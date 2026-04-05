import numpy as np
import pandas as pd

from settings import G
from src.models.solver import B2T1
from src.models.runge_kutta import rk4_step
from src.models.target_1 import f_2_center
from src.utils.utils import params_tuple
from src.utils.Exceptions import WrongEnergyError, WrongDistanceError
from itertools import product

A = 1/2
m_l, j_l, b_l = 2.1, 0.51, 0.36
m = 1
j = 1
b = 0.7

yC2 = [np.array([0, 0, 0]), np.array([10, 0, 0])]
params = params_tuple(m, j, b, A, yC2, 0)
y0 = np.array([5, 0, 0, 0., 0.322386, 0.0969581, 0, 0.48479, -1.44494])


t0_task = 0
tEnd_task = 50
n_task = 200
tau = (tEnd_task - t0_task) / n_task


solver = B2T1(f_2_center, rk4_step, t0_task, tEnd_task, tau, y0, params)
solver.solve()
solver.create_plotly_graph()
solver.create_graph_energy()
