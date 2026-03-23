import numpy as np

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

def Verle(f, t0, y0, tEnd, tau):
         def increment(f, t, y0, tau):
                  v = f(t,y0) - tau * y0
                  return tau * v
         t = []
         y = []
         t.append(t0)
         y.append(y0)


         while t0 < tEnd:
                  #print(y0) 
                  y0 = y0 + increment(f, t0, y0, tau) 
                  #print(increment(f, t0, y0, tau) )
                  #print(y0)
                  t0 = t0 + tau 
                  t.append(t0)
                  y.append(y0) 
             
         return np.array(t), np.array(y)

y0 =  np.array([5, 0, 0,    0., 0.322386, 0.0969581,    0, 0.48479, -1.44494])
 
m = 1
j = 1
b = 0.7
A = 1/2

yC2 = [0,0,0]

t0 = 0
tEnd = 100
n = tEnd * 1000
tau = (tEnd-t0)/n
t = np.linspace(t0,tEnd,n)

t, y = Verle(f, t0, y0, tEnd, tau)

print(y[-1])


# 3D ГРАФИКИ
import pandas as pd
from plotly.offline import plot
import plotly.graph_objs as go

fig = go.Figure()

# траектория таргета
df = pd.DataFrame({"x": y[:,0], "y": y[:,1], "z": y[:,2]})
df1 = df.drop(labels = [i for i in range(0,len(y[:,1])-1)],axis = 0)
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
                         size=30,
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
                         size=30,
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



# =============================================================================
# dfk1 = pd.DataFrame({"x": P1[:,0], "y": P1[:,1], "z": P1[:,2]})
# dfk11 = dfk1.drop(labels = [i for i in range(0,len(P1[:,1])-1)], axis = 0)
# 
# fig.add_trace(go.Scatter3d(x=dfk1['x'], y=dfk1['y'], z = dfk1['z'],
#                            mode = 'lines',
#                            name='траектория K1',
#                            line=dict(
#                                color=['black'],  
#                                colorscale='Viridis'),
#                            opacity=0.5))
# fig.add_trace(go.Scatter3d(x=dfk11['x'], y=dfk11['y'], z = dfk11['z'],
#                      mode='markers',
#                      name='Количество движения K1',
#                      marker=dict(
#                          size=8,
#                          color=['green'])))
# 
# 
# dfk2 = pd.DataFrame({"x": P2[:,0], "y": P2[:,1], "z": P2[:,2]})
# dfk22 = dfk2.drop(labels = [i for i in range(0,len(P1[:,1])-1)], axis = 0)
# 
# fig.add_trace(go.Scatter3d(x=dfk2['x'], y=dfk2['y'], z = dfk2['z'],
#                            mode = 'lines',
#                            name='траектория K2',
#                            line=dict(
#                                color=['black'],
#                                colorscale='Viridis'),
#                            opacity=0.5))
# fig.add_trace(go.Scatter3d(x=dfk22['x'], y=dfk22['y'], z = dfk22['z'],
#                      mode='markers',
#                      name='Собственный кинетический момент K2',
#                      marker=dict(
#                          size=8,
#                          color=['green'])))
# =============================================================================




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


plot(fig)
# fig.write_html('../../data/prob_v2.html')
