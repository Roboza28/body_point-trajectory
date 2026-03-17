import numpy as np


def f_1_center(y, params):
    m, j, b = params.m, params.j, params.b

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


def f_1_center_and_magnetic_field(y, params):
    m, j, b, A, B = params.m, params.j, params.b, params.A, params.B

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


def f_2_center(y, params):
    m, j, b, A, yC2 = params.m, params.j, params.b, params.A, params.yC2

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


def f_2_center_mat_point(y, params):
    pass

def f_2_center(y, params):
    m, j, b, A, yC2 = params.m, params.j, params.b, params.A, params.yC2
