import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
# from scipy.integrate import *
from matplotlib import pyplot as plt
from plotly.offline import plot
import plotly.graph_objs as go


def angle_between(v1, v2):
    dot_pr = v1.dot(v2)
    norms = np.linalg.norm(v1) * np.linalg.norm(v2)
    return np.rad2deg(np.arccos(dot_pr / norms))


def runge_kutta(f, t0, y0, tEnd, tau):
    def increment(f, t, y, tau):
        k1 = tau * f(t, y)
        k2 = tau * f(t + (1 / 4) * tau, y + (1 / 4) * k1)
        k3 = tau * f(t + (3 / 8) * tau, y + (3 / 32) * k1 + (9 / 32) * k2)
        k4 = tau * f(t + (12 / 13) * tau, y + (1932 / 2197) * k1 - (7200 / 2197) * k2 + (7296 / 2197) * k3)
        k5 = tau * f(t + tau, y + (439 / 216) * k1 - 8 * k2 + (3680 / 513) * k3 - (845 / 4104) * k4)
        k6 = tau * f(t + (1 / 2) * tau,
                     y - (8 / 27) * k1 + 2 * k2 - (3544 / 2565) * k3 + (1859 / 4104) * k4 - (11 / 40) * k5)
        return (16 / 135) * k1 + (6656 / 12825) * k3 + (28561 / 56430) * k4 - (9 / 50) * k5 + (2 / 55) * k6

    t = []
    y = []
    t.append(t0)
    y.append(y0)

    # E = []
    P1 = []
    P2 = []

    P3, P4 = [], []

    while t0 < tEnd:
        y0 = y0 + increment(f, t0, y0, tau)
        t0 = t0 + tau
        t.append(t0)
        y.append(y0)

        R = np.array([y0[0], y0[1], y0[2]])
        K1 = np.array([y0[3], y0[4], y0[5]])
        K2 = np.array([y0[6], y0[7], y0[8]])

        # e_theory = np.array([np.cos(gamma_theory), np.sin(gamma_theory), 0])
        # e_practice = np.array([np.cos(gamma_practice), np.sin(gamma_practice), 0])

        v = 1 / (m * j - b ** 2) * (j * K1 - b * K2)
        w = 1 / (b ** 2 - m * j) * (b * K1 - m * K2)

        # K = 1 / 2 * m * np.dot(v, v) + b * np.dot(v, w) + 1 / 2 * j * np.dot(w, w)
        # U = -A * (1 / (np.linalg.norm(R)))
        # E.append(K + U)

        K1N = K1 / np.linalg.norm(K1)
        P1.append(K1N)
        K2N = K2 / np.linalg.norm(K2)
        P2.append(K2N)
        # print(np.dot(R/np.linalg.norm(R), K2N))
        print(w)
        # P3.append(v / np.linalg.norm(v))
        P3.append(w / np.linalg.norm(w))
        # P4.append(np.dot(v, e_practice))

        # P11 = np.linalg.norm(w)
        # P22 = np.linalg.norm(v)
        # P1.append(P11)
        # P2.append(P22)

        print(np.linalg.norm(np.dot(K1, K1)))

    # return np.array(t), np.array(y), np.array(E), np.array(P1), np.array(P2)
    return np.array(t), np.array(y), np.array(P1), np.array(P2), np.array(P3) # , np.array(P4)


def f(t, y):
    f = np.zeros([9])
    f[0] = 1 / (m * j - b ** 2) * (j * y[3] - b * y[6])
    f[1] = 1 / (m * j - b ** 2) * (j * y[4] - b * y[7])
    f[2] = 1 / (m * j - b ** 2) * (j * y[5] - b * y[8])

    f[3] = - y[0] / ((y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** (3 / 2))
    f[4] = - y[1] / ((y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** (3 / 2))
    f[5] = - y[2] / ((y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** (3 / 2))

    f[6] = -b / (m * j - b ** 2) * (y[4] * y[8] - y[5] * y[7])
    f[7] = -b / (m * j - b ** 2) * (y[5] * y[6] - y[3] * y[8])
    f[8] = -b / (m * j - b ** 2) * (y[3] * y[7] - y[4] * y[6])

    return f


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


# y0 = np.array([5, 0, 0, 0., 0.322386, 0.0969581, 0, 0.48479, -1.44494])
# # y0 =  np.array([5, 0, 0,    0., 0.20, 0.0969581,    0, 0.48479, -1.44494])
# # y0 =  np.array([5, 0, 0,    0., 0., 0.,    0, 0., -1.61193])

# =============================================================================

# gamma_practice = math.atan(q * (m * j - b ** 2)/(b * k2y0**2)) + 2.85*math.pi/180
# gamma_practice = math.pi/2 - math.atan(q * (m * j - b ** 2)/(b * k2y0**2))
# gamma_practice = 0
# gamma_theory = math.atan(q * (m * j - b ** 2) / (b * k2y0 ** 2))

y0 = np.array([rx0, ry0, rz0, k1x0, k1y0, k1z0, k2x0, k2y0, k2z0])

t0 = 0
tEnd = 1
# tEnd = 500
n = tEnd * 100
tau = (tEnd - t0) / n
# t = np.linspace(t0, tEnd, n)

# t, y, E, P1, P2 = runge_kutta(f, t0, y0, tEnd, tau)
t, y, P1, P2, P3 = runge_kutta(f, t0, y0, tEnd, tau)


# fig = plt.figure()
# plt.title("Зависимость полной энергии от времени")
# # plt.xlabel('t')
# # plt.ylabel('E')
# plt.grid(True)
# plt.plot(t[1:len(t)], P4)
# plt.ylim([E[0]-0.1, E[0]+0.1]) # масштабирование графика
# plt.plot(t[1:len(t)], P1, 'k')
# plt.plot(t[1:len(t)], P2, 'k')
# plt.show()


# 3D ГРАФИКИ
df = pd.DataFrame({"x": y[:, 0], "y": y[:, 1], "z": y[:, 2]})
df1 = df.drop(labels=[i for i in range(0, len(y[:, 1]) - 1)], axis=0)

xC, yC, zC = [0], [0], [0]

df2 = pd.DataFrame({"x": xC, "y": yC, "z": zC})

fig = go.Figure()

# =============================================================================
# fig.add_trace(go.Scatter3d(x=df['x'], y=df['y'], z = df['z'],
#                            mode = 'lines',
#                            name='траектория',
#                            line=dict(
#                                color=['black'],  # set color to an array/list of desired values
#                                colorscale='Viridis',
#                                dash = 'dashdot'),
#                            opacity=0.5))
# 
# =============================================================================
fig.add_trace(go.Scatter3d(x=df['x'], y=df['y'], z=df['z'],
                           mode='lines',
                           name='траектория',
                           line=dict(
                               color=['black'],  # set color to an array/list of desired values
                               colorscale='Viridis'),
                           opacity=0.5))

fig.add_trace(go.Scatter3d(x=df1['x'], y=df1['y'], z=df1['z'],
                           mode='markers',
                           name='актуальное положение цели',
                           marker=dict(
                               size=12,
                               color=['red'])))

fig.add_trace(go.Scatter3d(x=df2['x'], y=df2['y'], z=df2['z'],
                           mode='markers',
                           name='массивный центр',
                           marker=dict(
                               size=30,
                               color=['yellow'],  # set color to an array/list of desired values
                               colorscale='Viridis')))

# =============================================================================
# dfk1 = pd.DataFrame({"x": P1[:, 0], "y": P1[:, 1], "z": P1[:, 2]})
#
# fig.add_trace(go.Scatter3d(x=dfk1['x'], y=dfk1['y'], z=dfk1['z'],
#                            mode='lines',
#                            name='траектория K1',
#                            line=dict(
#                                color=['black'],
#                                colorscale='Viridis'),
#                            opacity=0.5))
# fig.add_trace(go.Scatter3d(x=dfk1.iloc[-1:]['x'], y=dfk1.iloc[-1:]['y'], z=dfk1.iloc[-1:]['z'],
#                            mode='markers',
#                            name='Количество движения K1',
#                            marker=dict(
#                                size=8,
#                                color=['green'])))
# fig.add_trace(go.Scatter3d(x=dfk1.iloc[:1]['x'], y=dfk1.iloc[:1]['y'], z=dfk1.iloc[:1]['z'],
#                            mode='markers',
#                            name='Собственный кинетический момент K2',
#                            marker=dict(
#                                size=8,
#                                color=['darkgreen'])))


# dfk2 = pd.DataFrame({"x": P2[:, 0], "y": P2[:, 1], "z": P2[:, 2]})
# dfk22 = dfk2.drop(labels = [i for i in range(0,len(P1[:,1])-1)], axis = 0)
#
# fig.add_trace(go.Scatter3d(x=dfk2['x'], y=dfk2['y'], z=dfk2['z'],
#                            mode='lines',
#                            name='траектория K2',
#                            line=dict(
#                                color=['red'],
#                                colorscale='Viridis'),
#                            opacity=1))
# fig.add_trace(go.Scatter3d(x=dfk2.iloc[-1:]['x'], y=dfk2.iloc[-1:]['y'], z=dfk2.iloc[-1:]['z'],
#                            mode='markers',
#                            name='Собственный кинетический момент K2',
#                            marker=dict(
#                                size=8,
#                                color=['red'])))
# fig.add_trace(go.Scatter3d(x=dfk2.iloc[:1]['x'], y=dfk2.iloc[:1]['y'], z=dfk2.iloc[:1]['z'],
#                            mode='markers',
#                            name='Собственный кинетический момент K2',
#                            marker=dict(
#                                size=8,
#                                color=['darkred'])))


dfk3 = pd.DataFrame({"x": P3[:, 0], "y": P3[:, 1], "z": P3[:, 2]})

fig.add_trace(go.Scatter3d(x=dfk3['x'], y=dfk3['y'], z=dfk3['z'],
                           mode='lines',
                           name='траектория w',
                           line=dict(
                               color=['red'],
                               colorscale='Viridis'),
                           opacity=1))
fig.add_trace(go.Scatter3d(x=dfk3.iloc[-1:]['x'], y=dfk3.iloc[-1:]['y'], z=dfk3.iloc[-1:]['z'],
                           mode='markers',
                           name='Собственный кинетический момент w',
                           marker=dict(
                               size=8,
                               color=['red'])))
fig.add_trace(go.Scatter3d(x=dfk3.iloc[:1]['x'], y=dfk3.iloc[:1]['y'], z=dfk3.iloc[:1]['z'],
                           mode='markers',
                           name='Собственный кинетический момент w',
                           marker=dict(
                               size=8,
                               color=['darkred'])))

# =============================================================================

# print(gamma_theory*180/math.pi, gamma_practice*180/math.pi)
# length = 2

# fig.add_trace(go.Scatter3d(
#     x=[0, length * math.cos(gamma_theory)],  # Координаты по оси X
#     y=[0, length * math.sin(gamma_theory)],  # Координаты по оси Y
#     z=[0, 0],
#     mode='lines',
#     line=dict(
#            color=['black'],
#            colorscale='Viridis'),
#     opacity=1
# ))


# fig.add_trace(go.Scatter3d(
#     x=[0, length * math.cos(gamma_practice)],  # Координаты по оси X
#     y=[0, length * math.sin(gamma_practice)],  # Координаты по оси Y
#     z=[0, 0],
#     mode='lines',
#     line=dict(
#            color=['red'],
#            colorscale='Viridis'),
#     opacity=1
# ))

# plot(fig)
fig.write_html('../../data/trajectory_rk_1b_1t.html')
