import numpy as np
import matplotlib.pyplot as plt
#from scipy.integrate import *

def angle_between(v1, v2):
    dot_pr = v1.dot(v2)
    norms = np.linalg.norm(v1) * np.linalg.norm(v2)
    return np.rad2deg(np.arccos(dot_pr / norms))


def my_direvative(vector_act: np.array, vector_old: np.array, dt: float):
    
    vector_x = (vector_act[0] - vector_old[0]) / dt
    vector_y = (vector_act[1] - vector_old[1]) / dt
    vector_z = (vector_act[2] - vector_old[2]) / dt
    
    vector_new = np.array([vector_x, vector_y, vector_z])
    
    return vector_new

def rungeKutta(f, t0, y0, tEnd, tau):
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
                  
                  #if np.linalg.norm(R)>20:
                  #    print('Все плохо')
                  #    break
                  
                  v = 1/(m*j-b**2) * (j*K1 - b*K2)
                  w = 1/(b**2-m*j) * (b*K1 - m*K2)
                  
                  K = 1/2*m*np.dot(v,v) + b*np.dot(v,w) + 1/2*j*np.dot(w,w)
                  U = -A * (1/(np.linalg.norm(R)))\
                      -A * (1/(np.linalg.norm(Rs-R)))  
                  E.append(K+U)
                  
                  # R1 = Rs-R
                  #print(K,U)
                  # проверка величин на сохранность

                  K2q = np.cross(R,K1) + K2
                  #PK2q.append(np.array([K2q[0],K2q[1],K2q[2]]))

                  
                  #K2qP = np.dot(K2q,np.array([1,0,0]))
                  #P.append(K2qP)

                  #RRK1 = np.dot(np.cross(R,Rs),K1)
                  #RRK2 = np.dot(np.cross(R,Rs),K2)
                  #P.append(RRK1+RRK2)
                  
                  #RRw = (np.dot(np.cross(R,Rs),w))/(np.linalg.norm(Rs-R)**3)
                  #P.append(RRw)
                  #RRR=np.dot(np.cross(R,Rs),v)
                  #print(RRw)
                  #P.append(RRR)
                  
                  #af1 = (j/b)*np.dot(K1,Rs)
                  #af2 = -((m*j-b**2)/b)*np.dot(v,Rs)
                  #af3 = np.dot(np.cross(R,K1),Rs)
                  
                  #af11 = (j/b)*K1
                  #af12 = -((m*j-b**2)/b)*v
                  #af13 = np.cross(R,K1)
                  
                  #a = af11 + af12 + af13
                  
                  
                  f1 = -A*y0[0]/(((y0[0]**2+y0[1]**2+y0[2]**2))**(3/2)) - A*(y0[0]-yC2[0])/((((y0[0]-yC2[0])**2+(y0[1]-yC2[1])**2+(y0[2]-yC2[2])**2))**(3/2)) 
                  f2 = -A*y0[1]/(((y0[0]**2+y0[1]**2+y0[2]**2))**(3/2)) - A*(y0[1]-yC2[1])/((((y0[0]-yC2[0])**2+(y0[1]-yC2[1])**2+(y0[2]-yC2[2])**2))**(3/2)) 
                  f3 = -A*y0[2]/(((y0[0]**2+y0[1]**2+y0[2]**2))**(3/2)) - A*(y0[2]-yC2[2])/((((y0[0]-yC2[0])**2+(y0[1]-yC2[1])**2+(y0[2]-yC2[2])**2))**(3/2)) 
                  K1T = np.array([f1,f2,f3])
                  
                  f4 = -b/(m*j-b**2)*(y0[4]*y0[8] - y0[5]*y0[7])
                  f5 = -b/(m*j-b**2)*(y0[5]*y0[6] - y0[3]*y0[8])
                  f6 = -b/(m*j-b**2)*(y0[3]*y0[7] - y0[4]*y0[6])
                  K2T = np.array([f4,f5,f6])
                  K2qT = np.cross(R,K1T)
                  
                  f7 = 1/(m*j-b**2) * (j * y0[3] - b*y0[6])
                  f8 = 1/(m*j-b**2) * (j * y0[4] - b*y0[7])
                  f9 = 1/(m*j-b**2) * (j * y0[5] - b*y0[8])
                  RT = np.array([f7,f8,f9])
                  
                  R13 = 1/(np.linalg.norm(Rs-R)**3)
                  R3 = R13 + 1/(np.linalg.norm(R)**3) 
                  
                  
                  
                  #print(np.cross(v,R))
                  #print('*************************************')
                  #print(np.cross(K1T,R))
                  #print(np.cross(-A*R*R3+A*Rs*R13,R))
                  #print(np.cross(A*Rs*R13,R))
                  #print('*************************************')
                  
                  #print(np.dot(np.cross(K1T,R),Rs))
                  #print(np.dot(K2T,w))
                  
                  #print(np.dot(np.cross(RT,K1),v),np.dot(K2T,v),np.dot(K2qT,v))
                  #print(np.cross(RT,K1)+K2T,K2qT)
                  #print(np.dot(np.cross(R,K1T),v))
                  #print((np.dot(Rs,R))**2,np.dot(R,R))
                  
                  #a = np.cross(v,K1)+K2T
                  #b1 = K2qT
                  
                  #print(np.dot(np.cross(K1,RT),Rs),np.dot(np.cross(K1T,R),K1))
                  
                  #print(np.dot(np.cross(R,K1),Rs))
                  #print(np.dot(np.cross(R,v),Rs))
                  #print(w-(1/b*K1-m/b*v))
                  #print((np.dot(np.cross(R,Rs),w)))
                  
                  #print(np.dot(np.cross(K2q,R),Rs))
                  
                  #print(R)
                  #print(Rs)
                  #print(a)
                  #P.append(a)
                  #P.append(w/np.linalg.norm(w))
                  
                  #print((j/b)*K1 -((m*j-b**2)/b)*v + np.cross(R,K1))
                  

                  
                  #w1 = 1/(b**2-m*j) * (b*K1[0] - m*K2[0])
                  #w2 = 1/(b**2-m*j) * (b*K1[1] - m*K2[1])
                  #w3 = 1/(b**2-m*j) * (b*K1[2] - m*K2[2])
                  #WW = np.array([w1,w2,w3])
                  #P.append(WW)
                  
                  
                  # K1N = K1/np.linalg.norm(K1)
                  # P1.append(K1N)
                  # if np.linalg.norm(K2) != 0:
                  #     K2N = K2/np.linalg.norm(K2)
                  # else:
                  #     K2N = [0, 0, 0]
                  # P2.append(K2N)
                  
                  #R2 = np.dot(K2q,v)
                  #R2=np.dot(K2,Rs)
                  #P.append(R2)
                  
                  
                  #P1.append(np.linalg.norm(v))
                  #P2.append(np.linalg.norm(K1))
                  #P3.append(np.linalg.norm(K2))
                  
                  
                  #print(np.linalg.norm(K2q))
                  #print(angle_between(K1,K2))
                  #print(angle_between(K2q,Rs))
                  
                  
                  #print(np.linalg.norm(R)**3)
                  #print((R-Rs)^3)
                  
                  
                  # print((np.linalg.norm(R)**3)/(np.linalg.norm(R-Rs)**3))
                  
                  #print(np.dot(K2qT,R))
                  #print(np.dot(K2q,RT))
                  #print('-'*20)
                  if t0 > 2*tau:   
                      # v_xx = (R[0] - x_st) / tau
                      # v_yy = (R[1] - y_st) / tau
                      # v_zz = (R[2] - z_st) / tau
                      
                      # v_nnum = np.array([v_xx,v_yy,v_zz])
                      
                      print((np.linalg.norm(R) - np.linalg.norm(R_st)) / tau)#, np.linalg.norm(v))
                      # v_num = my_direvative(R, R_st, tau)
                      
                      # v1_num = my_direvative(Rs-R, R1_st, tau)
                      
                      # print(Rs-R,R1_st)
                      # print(np.linalg.norm(v_num),np.linalg.norm(v1_num))
                  
                  if t0>0:
                      # x_st = R[0]
                      # y_st = R[1]
                      # z_st = R[2]
                      R_st = R
                      R1_st = Rs-R
                    
                  
                  # print(f'K1: {np.linalg.norm(K1)}')
                  # print(f'K1t: {np.linalg.norm(K1T)}')
                  # print(f'R: {np.linalg.norm(R)}')
                  # print(f'R1: {np.linalg.norm(Rs-R)}')
                  # print(f'Rs: {np.linalg.norm(Rs)}')
                  # print(f'K2: {np.linalg.norm(K2)}')
                  # print(f'E: {np.linalg.norm(K+U)}')
                  
                  # print(f'v: {np.linalg.norm(v)}')
                  
                  # if t0 > 2*tau:  
                  #     print(f'v1: {np.linalg.norm(v1_num)}')
                  
                  print('-'*20)
                  #P1.append(np.linalg.norm(K1))
                  
         return np.array(t), np.array(y), np.array(E)# , np.array(P1), np.array(P2)# , np.array(P3)
     
        
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
     
y0 =  np.array([5, 0, 0,    0., 0.322386, 0.0969581,    0, 0.48479, -1.44494])
# y0 =  np.array([5, 0, 0,    0., 1, 0.0969581,    0, 0.48479, -1.44494])

# y0 =  np.array([1/2, 0, 0,    0, 1, 0,    0, 0, 0]) 





#y0 =  np.array([1/2, 0, -1,    0, 0.5 , 0.1,    0, 0, 0])  #!! отличное решение
#y0 =  np.array([1/2, 0, -1,    0, 0.1 , 0.1,    0, 0, 0])  #!! отличное решение
#y0 =  np.array([1/2, 0, 1,    0, 0.1 , 0.5,    0, 0, 0])  #!! отличное решение
#y0 =  np.array([1/2, 0, -1,    0, 0.01 , 0.1,    0, 0, 0])  #!! отличное решение



#y0 =  np.array([2/3, 0, 0,    0, 0.05 , 1,    0.1, 0, 0])  #!! отличное решение
    


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

t0 = 0
tEnd = 100
n = tEnd * 100
tau = (tEnd-t0)/n
t = np.linspace(t0,tEnd,n)


t, y, E = rungeKutta(f, t0, y0, tEnd, tau)
# t, y, E, P1, P2 = rungeKutta(f, t0, y0, tEnd, tau)
#print(E[0])
#print(E[-1])

fig = plt.figure()
plt.title("Зависимость полной энергии от времени")
plt.xlabel('t')
plt.ylabel('E')
plt.grid(True)
plt.axis([t[0],t[-1],(E[0]+E[-1])/2-0.1,(E[0]+E[-1])/2+0.1])
plt.plot(t[1:len(t)], E, 'k')
#plt.plot(t[1:len(t)], P, 'k')


#plt.plot(t[1:len(t)], P1, 'k')


#plt.plot(P1, P2, 'k')
#plt.plot(P1, P3, 'k')
#plt.plot(P2, P3, 'k')


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
