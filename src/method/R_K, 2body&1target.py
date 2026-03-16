import numpy as np
import matplotlib.pyplot as plt
#from scipy.integrate import *
import pandas as pd
# from plotly.offline import plot
import plotly.graph_objs as go

from src.method.runge_kutta import rk4_step
from src.models.solver import Solver
from src.models.target_1 import f_2_center
from src.utils.utils import params_tuple



y0 =  np.array([5, 0, 0,
                0., 0.322386, 0.0969581,
                0, 0.48479, -1.44494])
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


R1  = np.array([y0[0],y0[1],y0[2]])
K11 = np.array([y0[3],y0[4],y0[5]])
K21 = np.array([y0[6],y0[7],y0[8]])
#v1 = 1/(m*j-b**2) * (j*K11 - b*K21)
#print((j/b)*K11-((m*j-b**2)/b)*v1 + np.cross(R1,K11))

yC2 = [0.5,0,0]
#yC2 = [0,0,0]
#yC2 = [20,0,0]
#yC2 = [0.5,0,0]

t0_task = 0
tEnd_task = 100
n_task = 10000
tau = (tEnd_task - t0_task) / n_task


params = params_tuple(m, j, b, A, yC2)
t_sol, y_sol, E_sol, *other_data = Solver(f_2_center, rk4_step, t0_task, tEnd_task, tau, y0, params).solve()
# t_sol, y_sol, E_sol, *other_data = runge_kutta(f_2_center, t0_task, y0, tEnd_task, tau, params)
# t_sol, y_sol, E_sol, P1, P2 = runge_kutta(f, t0_task, y0, tEnd_task, tau)
# print(E[0], E[-1])

# fig_plt = plt.figure()
plt.title("Зависимость полной энергии от времени")
plt.xlabel('t')
plt.ylabel('E')
plt.grid(True)
plt.axis([t_sol[0], t_sol[-1], (E_sol[0] + E_sol[-1]) / 2 - 0.1, (E_sol[0] + E_sol[-1]) / 2 + 0.1])
plt.plot(t_sol[1:len(t_sol)], E_sol, 'k')
#plt.plot(t[1:len(t_sol)], P, 'k')
#plt.plot(t[1:len(t_sol)], P1, 'k')
#plt.plot(P1, P2, 'k')
#plt.plot(P1, P3, 'k')
#plt.plot(P2, P3, 'k')


# 3D ГРАФИКИ
fig = go.Figure()

# траектория тела-точки
df = pd.DataFrame({"x": y_sol[:,0], "y": y_sol[:,1], "z": y_sol[:,2]})
df1 = df.drop(labels = [i for i in range(0, len(y_sol[:,1]) - 1)], axis = 0)
fig.add_trace(go.Scatter3d(x=df['x'], y=df['y'], z = df['z'],
                           mode = 'lines',
                           name='траектория',
                           line=dict(
                               color=['black'],  # set color to an array/list of desired values
                               colorscale='Viridis'),
                           opacity=0.5))
fig.add_trace(go.Scatter3d(x=df1['x'], y=df1['y'], z = df1['z'],
                     mode='markers',
                     name='актуальное положение цели',
                     marker=dict(
                         size=15,
                         color=['red'])))



xC1 = [0]
yC1 = [0]
zC1 = [0]
df2 = pd.DataFrame({"x": xC1, "y": yC1, "z": zC1})
fig.add_trace(go.Scatter3d(x=df2['x'], y=df2['y'], z = df2['z'],
                     mode='markers',
                     name='массивный центр № 1',
                     marker=dict(
                         size=20,
                         color=['yellow'],  # set color to an array/list of desired values
                         colorscale='Viridis')))


xC3 = [yC2[0]]
yC3 = [yC2[1]]
zC3 = [yC2[2]]
df3 = pd.DataFrame({"x": xC3, "y": yC3, "z": zC3})
fig.add_trace(go.Scatter3d(x=df3['x'], y=df3['y'], z = df3['z'],
                     mode='markers',
                     name='массивный центр № 2',
                     marker=dict(
                         size=20,
                         color=['green'],  # set color to an array/list of desired values
                         colorscale='Viridis')))




# =============================================================================
# dfw1 = pd.DataFrame({"x": P[:,0], "y": P[:,1], "z": P[:,2]})
# #print(dfw1)
# dfw2 = dfw1.drop(labels = [i for i in range(0,len(P[:,1])-1)], axis = 0)
# print(dfw2)
# fig.add_trace(go.Scatter3d(x=dfw1['x'], y=dfw1['y'], z = dfw1['z'],
#                            mode = 'lines',
#                            name='траектория2',
#                            line=dict(
#                                color=['black'],  # set color to an array/list of desired values
#                                colorscale='Viridis'),
#                            opacity=0.5))
# fig.add_trace(go.Scatter3d(x=dfw2['x'], y=dfw2['y'], z = dfw2['z'],
#                      mode='markers',
#                      name='актуальное положение цели2',
#                      marker=dict(
#                          size=12,
#                          color=['blue'])))
# =============================================================================


'''
dfk1 = pd.DataFrame({"x": P1[:,0], "y": P1[:,1], "z": P1[:,2]})
dfk11 = dfk1.drop(labels = [i for i in range(0,len(P1[:,1])-1)], axis = 0)


fig.add_trace(go.Scatter3d(x=dfk1['x'], y=dfk1['y'], z = dfk1['z'],
                           mode = 'lines',
                           name='траектория K1',
                           line=dict(
                               color=['black'],  
                               colorscale='Viridis'),
                           opacity=0.5))
fig.add_trace(go.Scatter3d(x=dfk11['x'], y=dfk11['y'], z = dfk11['z'],
                     mode='markers',
                     name='Количество движения K1',
                     marker=dict(
                         size=8,
                         color=['purple'])))


dfk2 = pd.DataFrame({"x": P2[:,0], "y": P2[:,1], "z": P2[:,2]})
dfk22 = dfk2.drop(labels = [i for i in range(0,len(P2[:,1])-1)], axis = 0)

fig.add_trace(go.Scatter3d(x=dfk2['x'], y=dfk2['y'], z = dfk2['z'],
                           mode = 'lines',
                           name='траектория K2',
                           line=dict(
                               color=['black'],
                               colorscale='Viridis'),
                           opacity=0.5))
fig.add_trace(go.Scatter3d(x=dfk22['x'], y=dfk22['y'], z = dfk22['z'],
                     mode='markers',
                     name='Собственный кинетический момент K2',
                     marker=dict(
                         size=8,
                         color=['blue'])))
'''



# =============================================================================
# dfk2q = pd.DataFrame({"x": PK2q[:,0], "y": PK2q[:,1], "z": PK2q[:,2]})
# dfk22q = dfk2q.drop(labels = [i for i in range(0,len(PK2q[:,1])-1)], axis = 0)
# 
# fig.add_trace(go.Scatter3d(x=dfk2q['x'], y=dfk2q['y'], z = dfk2q['z'],
#                            mode = 'lines',
#                            name='траектория K2',
#                            line=dict(
#                                color=['black'],
#                                colorscale='Viridis'),
#                            opacity=0.5))
# fig.add_trace(go.Scatter3d(x=dfk22q['x'], y=dfk22q['y'], z = dfk22q['z'],
#                      mode='markers',
#                      name='Собственный кинетический момент K2',
#                      marker=dict(
#                          size=8,
#                          color=['green'])))
# =============================================================================

# plot(fig)
fig.write_html('../../data/trajectory_rk_2b_1t.html')
