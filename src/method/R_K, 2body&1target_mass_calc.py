import numpy as np
import pandas as pd
from multiprocessing import Pool
from settings import G
from src.models.solver import B2T1
from src.models.runge_kutta import rk4_step
from src.models.target_1 import f_2_center
from src.utils.utils import params_tuple
from src.utils.Exceptions import WrongEnergyError, WrongDistanceError
from itertools import product


# # Н.У.

# y0 =  np.array([5, 0, 0,    0., 0.322386, 0.0969581,   0, 0.48479, -1.44494])
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
yC2 = [np.array([0, 0, 0]), np.array([10, 0, 0])]


# yC2 = [np.array([0, 0, 0]), np.array([0.5, 0, 0])]
# yC2 = np.array([0.5, 0, 0])
#yC2 = [0,0,0]
#yC2 = [20,0,0]
#yC2 = [0.5,0,0]


# mc = np.sqrt(A / G)
# Rs = np.linalg.norm(yC2[1] - yC2[0])
# T = 2 * np.pi * np.sqrt(Rs ** 3 / (G * mc))
# print(T)

t0_task = 0
tEnd_task = 250
n_task = 200
tau = (tEnd_task - t0_task) / n_task


# params = params_tuple(m, j, b, A, yC2, 0)
# solver = B2T1(f_2_center, rk4_step, t0_task, tEnd_task, tau, y0, params)
# solver.solve()
# solver.create_plotly_graph()
# solver.create_graph_energy()




# list_with_m = np.concatenate([np.arange(0.01, 1.1, 0.05), np.arange(1.1, 10.1, 1)])
# list_with_j = np.concatenate([np.arange(0.01, 1.1, 0.05), np.arange(1.1, 10.1, 1)])
# list_with_b = np.concatenate([np.arange(0.01, 1.1, 0.05), np.arange(1.1, 10.1, 1)])
#
# combinations = list(product(list_with_m, list_with_j, list_with_b))
#
#
# result = []
# for m_l, j_l, b_l in combinations:
#     yC2 = [np.array([0, 0, 0]), np.array([10, 0, 0])]
#     params = params_tuple(m_l, j_l, b_l, A, yC2, 0)
#     y0 =  np.array([5, 0, 0,    0., 0.322386, 0.0969581,   0, 0.48479, -1.44494])
#
#     if np.allclose(m_l * j_l, b_l ** 2, rtol=0.01, atol=0.0):
#         continue
#
#     solver = B2T1(f_2_center, rk4_step, t0_task, tEnd_task, tau, y0, params)
#
#     try:
#         solver.solve()
#
#     except WrongEnergyError:
#         # result.append([m_l, j_l, b_l, 0])
#         continue
#
#     except WrongDistanceError:
#         # result.append([m_l, j_l, b_l, 1])
#         continue
#
#     except ZeroDivisionError:
#         # result.append([m_l, j_l, b_l, 2])
#         continue
#
#     # result.append([m_l, j_l, b_l])
#
#
# df = pd.DataFrame(result, columns=['m', 'j', 'b', 'tag'])


r0y  = np.arange(0, 11, 1)
k10x = np.concatenate([np.arange(0, 1, 0.1), np.arange(1, 11, 1)])
k10y = np.concatenate([np.arange(0, 1, 0.1), np.arange(1, 11, 1)])
k10z = np.concatenate([np.arange(0, 1, 0.1), np.arange(1, 11, 1)])
k20x = np.concatenate([np.arange(0, 1, 0.1), np.arange(1, 11, 1)])
k20y = np.concatenate([np.arange(0, 1, 0.1), np.arange(1, 11, 1)])
k20z = np.concatenate([np.arange(0, 1, 0.1), np.arange(1, 11, 1)])
# r0y  = np.arange(1, 10, 1)
# k10x = np.arange(0, 1, 0.1)
# k10y = np.arange(0, 1, 0.1)
# k10z = np.arange(0, 1, 0.1)
# k20x = np.arange(0, 1, 0.1)
# k20y = np.arange(0, 1, 0.1)
# k20z = np.arange(0, 1, 0.1)
combinations = list(product(r0y, k10x, k10y, k10z, k20x, k20y, k20z))
# combinations = list(product(r0y, k10y, k10z, k20y, k20z))


# for r0y_l, k10x_l, k10y_l, k10z_l, k20x_l, k20y_l, k20z_l in combinations:
#     y0 =  np.array([5, r0y_l, 0,    k10x_l, k10y_l, k10z_l,   k20x_l, k20y_l, k20z_l])
#
#     solver = B2T1(f_2_center, rk4_step, t0_task, tEnd_task, tau, y0, params)
#
#     try:
#         solver.solve()
#
#     except WrongEnergyError:
#         # result.append([m_l, j_l, b_l, 0])
#         continue
#
#     except WrongDistanceError:
#         # result.append([m_l, j_l, b_l, 1])
#         continue
#
#     except ZeroDivisionError:
#         # result.append([m_l, j_l, b_l, 2])
#         continue
#
#     result.append([r0y_l, k10x_l, k10y_l, k10z_l, k20x_l, k20y_l, k20z_l])
#
#     i+=1
#
#     print(i / e * 100)


params = params_tuple(m, j, b, A, yC2, 0)
def process_chunk(chunk_data):
    r0y_l, k10x_l, k10y_l, k10z_l, k20x_l, k20y_l, k20z_l = chunk_data
    # r0y_l, k10y_l, k10z_l, k20y_l, k20z_l = chunk_data
    y0 =  np.array([5, r0y_l, 0,    0, k10y_l, k10z_l,   0, k20y_l, k20z_l])
    solver = B2T1(f_2_center, rk4_step, t0_task, tEnd_task, tau, y0, params)

    try:
        solver.solve()

    except WrongEnergyError:
        # result.append([m_l, j_l, b_l, 0])
        # return 0
        return r0y_l, k10x_l, k10y_l, k10z_l, k20x_l, k20y_l, k20z_l, 0
        # return r0y_l, k10y_l, k10z_l, k20y_l, k20z_l, 0

    except WrongDistanceError:
        # result.append([m_l, j_l, b_l, 1])
        # return 0
        return r0y_l, k10x_l, k10y_l, k10z_l, k20x_l, k20y_l, k20z_l, 1
        # return r0y_l, k10y_l, k10z_l, k20y_l, k20z_l, 1

    except ZeroDivisionError:
        # result.append([m_l, j_l, b_l, 2])
        # return 0
        return r0y_l, k10x_l, k10y_l, k10z_l, k20x_l, k20y_l, k20z_l, 2
        # return r0y_l, k10y_l, k10z_l, k20y_l, k20z_l, 2

    return r0y_l, k10x_l, k10y_l, k10z_l, k20x_l, k20y_l, k20z_l, 3
    # return r0y_l, k10y_l, k10z_l, k20y_l, k20z_l, 3


if __name__ == "__main__":
    results = []
    print(0)
    with Pool(8) as pool:
        # pool.map автоматически распределит элементы по 4 процессам
        batch_results = pool.map(process_chunk, combinations)
        results.extend(batch_results)

    # df = pd.DataFrame(results, columns=['r0', 'k1x', 'k1y', 'k1z', 'k2x', 'k2y', 'k2z', 'tag'])
    df = pd.DataFrame(results, columns=['r0', 'k1y', 'k1z', 'k2y', 'k2z', 'tag'])
    df = df[df['tag'] == 3]

    df.to_excel('../../data/result5.xlsx', index=False)
