import numpy as np
from scipy.integrate import odeint



def model(y,t):
    xRdt = 1/(m*j-b**2) * (j * y[3] - b*y[6])
    yRdt = 1/(m*j-b**2) * (j * y[4] - b*y[7])
    zRdt = 1/(m*j-b**2) * (j * y[5] - b*y[8])
    
    xK1dt = -A*y[0]/(((y[0]**2+y[1]**2+y[2]**2))**(3/2)) 
    yK1dt = -A*y[1]/(((y[0]**2+y[1]**2+y[2]**2))**(3/2)) 
    zK1dt = -A*y[2]/(((y[0]**2+y[1]**2+y[2]**2))**(3/2)) 
    
    xK2dt = -b/(m*j-b**2)*(y[4]*y[8] - y[5]*y[7])
    yK2dt = -b/(m*j-b**2)*(y[5]*y[6] - y[3]*y[8])
    zK2dt = -b/(m*j-b**2)*(y[3]*y[7] - y[4]*y[6])
    
    return [xRdt,yRdt,zRdt,xK1dt,yK1dt,zK1dt,xK2dt,yK2dt,zK2dt]

y0 = [5, 0, 0,    0., 0.322386, 0.0969581,    0, 0.48479, -1.44494]

m = 1
j = 1
b = 0.7
A = 1

t0 = 0
tEnd = 10
n = 100
tau = (tEnd-t0)/n
t = np.linspace(t0,tEnd,n)


#sol = odeint(model,y0,t)
# sol[:,1] массив со скоростями (z)
# sol[:,0] массив с координатами (y)

x = []
y1 = []
z = []
E = []


for i in range(1,n):
    tspan = [t[i-1],t[i]]
    y = odeint(model,y0,tspan)
    x.append(y[1][0])
    y1.append(y[1][1])
    z.append(y[1][2])
    
    
    K1 = np.array([y[1][3],y[1][4],y[1][5]])
    K2 = np.array([y[1][6],y[1][7],y[1][8]])

    v = 1/(m*j-b**2) * (j*K1 - b*K2)
    w = 1/(b**2-m*j) * (b*K1 - m*K2)
    
    K = 1/2*m*np.dot(v,v) + b*np.dot(v,w) + 1/2*j*np.dot(w,w)
    U = -A*(1/(np.linalg.norm(np.array([y[1][0],y[1][1],y[1][2]]))))
    
    E.append(K+U)
    
    
    y0 = y[1]


import matplotlib.pyplot as plt
fig = plt.figure()
plt.title("Зависимость полной энергии от времени")
plt.xlabel('t')
plt.ylabel('E')
plt.grid(True)
plt.plot(t[1:len(t)], E, 'k')

import pandas as pd
from plotly.offline import plot
import plotly.graph_objs as go



df = pd.DataFrame({"x": x, "y": y1, "z": z})
print(df)
df1 = df.drop(labels = [i for i in range(0,len(y1)-1)],axis = 0)


xC = [0]
yC = [0]
zC = [0]

df2 = pd.DataFrame({"x": xC, "y": yC, "z": zC})


fig = go.Figure()
fig.add_trace(go.Scatter3d(x=df['x'], y=df['y'], z = df['z'],
                           mode = 'lines',
                           name='траектория',
                           marker=dict(
                               color=['black'])))

fig.add_trace(go.Scatter3d(x=df1['x'], y=df1['y'], z = df1['z'],
                     mode='markers',
                     name='актуальное положение цели',
                     marker=dict(
                         size=12,
                         color=['red'],  # set color to an array/list of desired values
                         colorscale='Viridis')))


fig.add_trace(go.Scatter3d(x=df2['x'], y=df2['y'], z = df2['z'],
                     mode='markers',
                     name='массивный центр',
                     marker=dict(
                         size=30,
                         color=['yellow'],  # set color to an array/list of desired values
                         colorscale='Viridis')))


# plot(fig)
fig.write_html('../../data/trajectory_odeint_1b_1t.html')
