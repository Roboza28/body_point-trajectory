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
    m, j, b, A = params.m, params.j, params.b, params.A
    Rs1, Rs2 = params.coord_centers
    # yC2 = params.coord_centers

    K1 = np.array([y[3], y[4], y[5]])
    K2 = np.array([y[6], y[7], y[8]])
    R  = np.array([y[0], y[1], y[2]])

    s = np.zeros([9])
    s[0], s[1], s[2] = 1/(m*j-b**2) * (j * K1 - b * K2)
    s[3], s[4], s[5] = - A / np.linalg.norm(R - Rs1)**3 * (R - Rs1)  - A / np.linalg.norm(R - Rs2)**3 * (R - Rs2)
    s[6], s[7], s[8] = -b / (m * j - b ** 2) * np.cross(K1, K2)

    return s


def f_2_center_mat_point(y, params):
    m, j, b, A = params.m, params.j, params.b, params.A
    Rs1, Rs2 = params.coord_centers

    R = np.array([y[0], y[2], y[4]])
    v = np.array([y[1], y[3], y[5]])

    s = np.zeros([6])
    s[0], s[2], s[4] =  v
    s[1], s[3], s[5] = - A / np.linalg.norm(R - Rs1)**3 * (R - Rs1)  - A / np.linalg.norm(R - Rs2)**3 * (R - Rs2)

    f = np.zeros([6])
    f[0] = y[1]
    f[1] = -A * (y[0]) / ((y[0] ** 2 + y[2] ** 2 + y[4] ** 2) ** (3 / 2)) \
           - A * (y[0] - Rs2[0]) / (((y[0] - Rs2[0]) ** 2 + (y[2] - Rs2[1]) ** 2 + (y[4] - Rs2[2]) ** 2) ** (3 / 2))
    f[2] = y[3]
    f[3] = -A * (y[2]) / ((y[0] ** 2 + y[2] ** 2 + y[4] ** 2) ** (3 / 2)) \
           - A * (y[2] - Rs2[1]) / (((y[0] - Rs2[0]) ** 2 + (y[2] - Rs2[1]) ** 2 + (y[4] - Rs2[2]) ** 2) ** (3 / 2))
    f[4] = y[5]
    f[5] = -A * (y[4]) / ((y[0] ** 2 + y[2] ** 2 + y[4] ** 2) ** (3 / 2)) \
           - A * (y[4] - Rs2[2]) / (((y[0] - Rs2[0]) ** 2 + (y[2] - Rs2[1]) ** 2 + (y[4] - Rs2[2]) ** 2) ** (3 / 2))

    return s
