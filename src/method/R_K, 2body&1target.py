import numpy as np
import matplotlib.pyplot as plt
#from scipy.integrate import *
import pandas as pd
# from plotly.offline import plot
import plotly.graph_objs as go


def angle_between(v1, v2):
    dot_pr = v1.dot(v2)
    norms = np.linalg.norm(v1) * np.linalg.norm(v2)
    return np.rad2deg(np.arccos(dot_pr / norms))


def derivative_vector(vector_act: np.ndarray, vector_old: np.ndarray, dt: float):
    vector_x = (vector_act[0] - vector_old[0]) / dt
    vector_y = (vector_act[1] - vector_old[1]) / dt
    vector_z = (vector_act[2] - vector_old[2]) / dt
    
    vector_new = np.array([vector_x, vector_y, vector_z])

    return vector_new


def runge_kutta(f, t0, y0, tEnd, tau):
    def increment(f, t, y, tau):
          k1 = tau * f(t,y)
          k2 = tau * f(t+(1/4)*tau,y+(1/4)*k1)
          k3 = tau * f(t+(3/8)*tau,y+(3/32)*k1+(9/32)*k2)
          k4 = tau * f(t+(12/13)*tau,y+(1932/2197)*k1-(7200/2197)*k2+(7296/2197)*k3)
          k5 = tau * f(t+tau,y+(439/216)*k1-8*k2+(3680/513)*k3 -(845/4104)*k4)
          k6 = tau * f(t+(1/2)*tau,y-(8/27)*k1+2*k2-(3544/2565)*k3 +(1859/4104)*k4-(11/40)*k5)
          return (16/135)*k1+(6656/12825)*k3+(28561/56430)*k4-(9/50)*k5+(2/55)*k6


    t = []
    y = []
    t.append(t0)
    y.append(y0)

    E = []
    # P = []
    # P1 = []
    # P2 = []
    # P3 = []
    #PK2q = []

    while t0 < tEnd:

        y0 = y0 + increment(f, t0, y0, tau)
        t0 = t0 + tau
        t.append(t0)
        y.append(y0)

        Rs = np.array([yC2[0],yC2[1],yC2[2]])
        R  = np.array([y0[0],y0[1],y0[2]])
        K1 = np.array([y0[3],y0[4],y0[5]])
        K2 = np.array([y0[6],y0[7],y0[8]])

        K1t = -A / np.linalg.norm(R) ** 3 * R - A / np.linalg.norm(R - Rs) ** 3 * (R - Rs)
        K2t = -b / (m * j - b ** 2) * np.cross(K1, K2)
        Rt  =  1 / (m * j - b ** 2) * (j * K1 - b * K2)

        v = 1/(m*j-b**2) * (j * K1 - b * K2)
        w = 1/(b**2-m*j) * (b * K1 - m * K2)

        K = 1/2 * m * np.dot(v,v) + b * np.dot(v,w) + 1/2 * j * np.dot(w,w)
        U = -A / np.linalg.norm(R) - A / np.linalg.norm(R - Rs)
        E.append(K + U)
        # print(K + U)

        K2q = np.cross(R,K1) + K2
        R13 = 1/(np.linalg.norm(Rs-R)**3)
        R3 = R13 + 1/(np.linalg.norm(R)**3)


        # K1N = K1/np.linalg.norm(K1)
        # P1.append(K1N)
        # if np.linalg.norm(K2) != 0:
        #     K2N = K2/np.linalg.norm(K2)
        # else:
        #     K2N = [0, 0, 0]
        # P2.append(K2N)

        #P1.append(np.linalg.norm(v))
        #P2.append(np.linalg.norm(K1))
        #P3.append(np.linalg.norm(K2))


        # print(np.dot(K1N, K2N))
        # print(np.dot(np.cross(R, K1), v + b/(m*j-b**2)*K2))
        # print(np.linalg.norm(v + b/(m*j-b**2)*K2))
        # print(np.dot(v + b/(m*j-b**2)*K2, v + b/(m*j-b**2)*K2))
        # print(np.dot(np.cross(R, K1), v + b/(m*j-b**2)*K2))
        # print(np.dot(v, v)+2*b*j/(m*j-b**2)**2*np.dot(K1, K2))
        # print(np.dot(K2t, K1))
        # print(np.dot(K1T, v+np.cross(w, R)))
        # print(np.dot(np.cross(R, K1T), w) + np.dot(K1T, v))

        # print('-'*20)
        #P1.append(np.linalg.norm(K1))

    return np.array(t), np.array(y), np.array(E) # , np.array(P1), np.array(P2)# , np.array(P3)


def f(t, y):
    f = np.zeros([9])

    f[0] = 1/(m*j-b**2) * (j * y[3] - b*y[6])
    f[1] = 1/(m*j-b**2) * (j * y[4] - b*y[7])
    f[2] = 1/(m*j-b**2) * (j * y[5] - b*y[8])

    f[3] = -A*y[0]/(((y[0]**2+y[1]**2+y[2]**2))**(3/2)) - A*(y[0]-yC2[0])/((((y[0]-yC2[0])**2+(y[1]-yC2[1])**2+(y[2]-yC2[2])**2))**(3/2))
    f[4] = -A*y[1]/(((y[0]**2+y[1]**2+y[2]**2))**(3/2)) - A*(y[1]-yC2[1])/((((y[0]-yC2[0])**2+(y[1]-yC2[1])**2+(y[2]-yC2[2])**2))**(3/2))
    f[5] = -A*y[2]/(((y[0]**2+y[1]**2+y[2]**2))**(3/2)) - A*(y[2]-yC2[2])/((((y[0]-yC2[0])**2+(y[1]-yC2[1])**2+(y[2]-yC2[2])**2))**(3/2))

    f[6] = -b/(m*j-b**2)*(y[4]*y[8] - y[5]*y[7])
    f[7] = -b/(m*j-b**2)*(y[5]*y[6] - y[3]*y[8])
    f[8] = -b/(m*j-b**2)*(y[3]*y[7] - y[4]*y[6])

    return f


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

yC2 = [0.1,0,0]
#yC2 = [0,0,0]
#yC2 = [20,0,0]
#yC2 = [0.5,0,0]

t0_task = 0
tEnd_task = 100
n_task = tEnd_task * 100
tau = (tEnd_task - t0_task) / n_task
t = np.linspace(t0_task, tEnd_task, n_task)


t_sol, y_sol, E_sol = runge_kutta(f, t0_task, y0, tEnd_task, tau)
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
