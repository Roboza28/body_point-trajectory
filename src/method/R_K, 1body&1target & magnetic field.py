import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import *


def rungeKutta(f, t0, y0, tEnd, tau):
         ii = 0
         def increment(f, t, y, tau):
                  k1 = tau  *f(t,y)
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

         while t0 < tEnd:
                  y0 = y0 + increment(f, t0, y0, tau) 
                  t0 = t0 + tau 
                  t.append(t0)
                  y.append(y0) 
                  
                  K1 = np.array([y0[3],y0[4],y0[5]])
                  K2 = np.array([y0[6],y0[7],y0[8]])

                  v = 1/(m*j-b**2) * (j*K1 - b*K2)
                  w = 1/(b**2-m*j) * (b*K1 - m*K2)
    
                  K = 1/2*m*np.dot(v,v) + b*np.dot(v,w) + 1/2*j*np.dot(w,w)
                  
                  U = -A*(1/(np.linalg.norm(np.array([y0[0],y0[1],y0[2]]))))

                  E.append(K+U)
                  

         return np.array(t), np.array(y), np.array(E)
     
        
def f(t, y):    
         f = np.zeros([9])
         f[0] = 1/(m*j-b**2) * (j * y[3] - b*y[6])
         f[1] = 1/(m*j-b**2) * (j * y[4] - b*y[7])
         f[2] = 1/(m*j-b**2) * (j * y[5] - b*y[8])
    
         # с постоянным магнитным полем, направленным по оси Z
         f[3] = -A*y[0]/(((y[0]**2+y[1]**2+y[2]**2))**(3/2)) + A*f[1]*B
         f[4] = -A*y[1]/(((y[0]**2+y[1]**2+y[2]**2))**(3/2)) - A*f[0]*B
         f[5] = -A*y[2]/(((y[0]**2+y[1]**2+y[2]**2))**(3/2)) 
         

         f[6] = -b/(m*j-b**2)*(y[4]*y[8] - y[5]*y[7])
         f[7] = -b/(m*j-b**2)*(y[5]*y[6] - y[3]*y[8])
         f[8] = -b/(m*j-b**2)*(y[3]*y[7] - y[4]*y[6])
         return f
     
y0 =  np.array([5, 0, 0,    0., 0.322386, 0.0969581,    0, 0.48479, -1.44494])


m = 1
j = 1
b = 0.7
A = 1
B = 0.1

t0 = 0
tEnd = 500
n = tEnd*100
tau = (tEnd-t0)/n
t = np.linspace(t0,tEnd,n)

t, y, E = rungeKutta(f, t0, y0, tEnd, tau)

fig = plt.figure()
plt.title("Зависимость полной энергии от времени")
plt.xlabel('t')
plt.ylabel('E')
plt.grid(True)
plt.plot(t[1:len(t)], E, 'k')
plt.ylim([2.16, 2.17]) # масштабирование графика


# 3D ГРАФИКИ
import pandas as pd
from plotly.offline import plot
import plotly.graph_objs as go



df = pd.DataFrame({"x": y[:,0], "y": y[:,1], "z": y[:,2]})
df1 = df.drop(labels = [i for i in range(0,len(y[:,1])-1)],axis = 0)

xC = [0]
yC = [0]
zC = [0]

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
                         size=12,
                         color=['red'])))

fig.add_trace(go.Scatter3d(x=df2['x'], y=df2['y'], z = df2['z'],
                     mode='markers',
                     name='массивный центр',
                     marker=dict(
                         size=30,
                         color=['yellow'],  # set color to an array/list of desired values
                         colorscale='Viridis')))

# plot(fig)
fig.write_html('../../data/trajectory_rk_1b_1t_magnetic.html')
