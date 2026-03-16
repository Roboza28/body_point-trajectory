import numpy as np
import pandas as pd


a, b = 2, 1

df = pd.read_excel('..\..\data\Ellipse.xlsx')
x, y = df.x, df.y

c = np.sqrt(a**2 - b**2)
print(c)

F1 = np.array([ c, 0])
F2 = np.array([-c, 0])

for i, j in zip(x, y):
    act_point = np.array([i, j])
    R0 = np.linalg.norm(act_point - F1)
    R1 = np.linalg.norm(act_point - F2)
    # print(R0 + R1)
    # print(R0 ** 3 + R1 ** 3)
    print(R0, R1)