import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import *
import pandas as pd
from plotly.offline import plot
import plotly.graph_objs as go


def rungeKutta(f, t0, y0, tEnd, tau):
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
    E = []
    t.append(t0)
    y.append(y0)

    while t0 < tEnd:
        y0 = y0 + increment(f, t0, y0, tau)
        t0 = t0 + tau
        t.append(t0)
        y.append(y0)

        Rs = np.array([yC2[0], yC2[1], yC2[2]])
        # R = np.array([y0[1], y0[3], y0[5]])
        # v = np.array([y0[0], y0[2], y0[4]])

        R = np.array([y0[0], y0[2], y0[4]])
        v = np.array([y0[1], y0[3], y0[5]])

        # if np.linalg.norm(R)>20:
        #    print('Все плохо')
        #    break

        K = 1 / 2 * m * np.dot(v, v)
        U = -A * (1 / (np.linalg.norm(R)) + 1 / (np.linalg.norm(R - Rs)))
        # print('-'*10)
        # print(K + U)
        print('-'*10)
        # E.append(K + U)
        r = np.linalg.norm(R)
        r1 = np.linalg.norm(R - Rs)
        S = 5*(1/r+1/r1) + 1/r1**3 + 1/r**3 - r**2/r1**3 - r1**2/r**3
        print(S)
        E.append(S)
    return np.array(t), np.array(y)# , np.array(E)


def f(t, y):
    f = np.zeros([6])
    f[0] = y[1]
    f[1] = -A * (y[0]) / ((y[0] ** 2 + y[2] ** 2 + y[4] ** 2) ** (3 / 2)) \
           - A * (y[0] - yC2[0]) / (((y[0] - yC2[0]) ** 2 + (y[2] - yC2[1]) ** 2 + (y[4] - yC2[2]) ** 2) ** (3 / 2))
    f[2] = y[3]
    f[3] = -A * (y[2]) / ((y[0] ** 2 + y[2] ** 2 + y[4] ** 2) ** (3 / 2)) \
           - A * (y[2] - yC2[1]) / (((y[0] - yC2[0]) ** 2 + (y[2] - yC2[1]) ** 2 + (y[4] - yC2[2]) ** 2) ** (3 / 2))
    f[4] = y[5]
    f[5] = -A * (y[4]) / ((y[0] ** 2 + y[2] ** 2 + y[4] ** 2) ** (3 / 2)) \
           - A * (y[4] - yC2[2]) / (((y[0] - yC2[0]) ** 2 + (y[2] - yC2[1]) ** 2 + (y[4] - yC2[2]) ** 2) ** (3 / 2))

    return f


m = 1

A = 1/2

# X0 = 2, VY0 = 1.5

x0 = 1
y0 = 0
z0 = 0
vx0 = 0
vy0 = 1.1
vz0 = 0

y0 = np.array([x0, vx0, y0, vy0, z0, vz0])
yC2 = [0.1, 0, 0]

t0 = 0
tEnd = 1000
n = tEnd * 100
tau = (tEnd - t0) / n

t, y, E = rungeKutta(f, t0, y0, tEnd, tau)

fig = plt.figure()
plt.title("Зависимость полной энергии от времени")
plt.xlabel('t')
plt.ylabel('E')
plt.grid(True)
plt.plot(t[1:len(t)], E, 'k')
plt.show()

# plt.axis([t[0],t[-1],(E[0]+E[-1])/2-0.1,(E[0]+E[-1])/2+0.1])
plt.plot(t[1:len(t)], E, 'k')

plt.plot(t[0:len(t)], y[:, 0], 'k')

# 3D ГРАФИКИ

fig = go.Figure()

df = pd.DataFrame({"x": y[:, 0], "y": y[:, 2], "z": y[:, 4]})
dff = pd.DataFrame({"t": [i for i in range(0, len(y[:, 1]))]})
df[['t']] = dff
df1 = df.drop(labels=[i for i in range(0, len(y[:, 1]) - 1)], axis=0)
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

xC1 = [0]
yC1 = [0]
zC1 = [0]
df2 = pd.DataFrame({"x": xC1, "y": yC1, "z": zC1})
fig.add_trace(go.Scatter3d(x=df2['x'], y=df2['y'], z=df2['z'],
                           mode='markers',
                           name='массивный центр № 1',
                           marker=dict(
                               size=30,
                               color=['yellow'],  # set color to an array/list of desired values
                               colorscale='Viridis')))

xC3 = [yC2[0]]
yC3 = [yC2[1]]
zC3 = [yC2[2]]
df3 = pd.DataFrame({"x": xC3, "y": yC3, "z": zC3})
fig.add_trace(go.Scatter3d(x=df3['x'], y=df3['y'], z=df3['z'],
                           mode='markers',
                           name='массивный центр № 2',
                           marker=dict(
                               size=30,
                               color=['green'],  # set color to an array/list of desired values
                               colorscale='Viridis')))

# plot(fig)
fig.write_html('../../data/trajectory_mat_point.html')
