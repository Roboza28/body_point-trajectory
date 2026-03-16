import numpy as np


class Solver:
    def __init__(self, system_of_task, inkrement, t_start, t_end, t_delta, y0, params_task):
        self.system_of_task = system_of_task
        self.inkrement = inkrement
        self.t_start = t_start
        self.t_end = t_end
        self.t_delta = t_delta
        self.params_task = params_task
        self.y0 = y0

    def solve(self):
        t = []
        y = []
        t.append(self.t_start)
        y.append(self.y0)

        E = []

        t_n = self.t_start
        y_n = self.y0

        while t_n < self.t_end:
            y_n1 = y_n + self.inkrement(self.system_of_task, y_n, self.t_delta, self.params_task)
            t_n1 = t_n + self.t_delta

            y.append(y_n1)
            t.append(t_n1)


            Rs = np.array([self.params_task.yC2[0], self.params_task.yC2[1], self.params_task.yC2[2]])
            R  = np.array([y_n1[0], y_n1[1], y_n1[2]])
            K1 = np.array([y_n1[3], y_n1[4], y_n1[5]])
            K2 = np.array([y_n1[6], y_n1[7], y_n1[8]])

            K1t = -self.params_task.A / np.linalg.norm(R) ** 3 * R - self.params_task.A / np.linalg.norm(R - Rs) ** 3 * (R - Rs)
            K2t = -self.params_task.b / (self.params_task.m * self.params_task.j - self.params_task.b ** 2) * np.cross(K1, K2)
            Rt = 1 / (self.params_task.m * self.params_task.j - self.params_task.b ** 2) * (self.params_task.j * K1 - self.params_task.b * K2)

            # K1t = - q / np.linalg.norm(R)**3 * R - q / np.linalg.norm(R - Rs)**3
            # K2t = -b / (m * j - b ** 2) * np.cross(K1, K2)

            # e_theory = np.array([np.cos(gamma_theory), np.sin(gamma_theory), 0])
            # e_practice = np.array([np.cos(gamma_practice), np.sin(gamma_practice), 0])

            v = 1 / (self.params_task.m * self.params_task.j - self.params_task.b ** 2) * (self.params_task.j * K1 - self.params_task.b * K2)
            w = 1 / (self.params_task.b ** 2 - self.params_task.m * self.params_task.j) * (self.params_task.b * K1 - self.params_task.m * K2)


            K = 1 / 2 * self.params_task.m * np.dot(v, v) + self.params_task.b * np.dot(v, w) + 1 / 2 * self.params_task.j * np.dot(w, w)
            # U = -self.params_task.A * (1 / (np.linalg.norm(R)))
            U = -self.params_task.A / np.linalg.norm(R) - self.params_task.A / np.linalg.norm(R - Rs)
            E.append(K + U)

            # K1N = K1 / np.linalg.norm(K1)
            # P1.append(K1N)
            # K2N = K2 / np.linalg.norm(K2)
            # P2.append(K2N)

            # print(np.dot(K1N, K2N))
            # print(np.dot(np.cross(R, K1), v + b/(m*j-b**2)*K2))
            # print(np.linalg.norm(v + b/(m*j-b**2)*K2))
            # print(np.dot(v + b/(m*j-b**2)*K2, v + b/(m*j-b**2)*K2))
            # print(np.dot(v, v)+2*b*j/(m*j-b**2)**2*np.dot(K1, K2))
            # print(np.dot(K2t, K1))

            # print(np.dot(R/np.linalg.norm(R), K2N))
            # print(w)
            # P3.append(v / np.linalg.norm(v))
            # P3.append(w / np.linalg.norm(w))
            # P4.append(np.dot(v, e_practice))

            # P11 = np.linalg.norm(w)
            # P22 = np.linalg.norm(v)
            # P1.append(P11)
            # P2.append(P22)

            # print(np.linalg.norm(np.dot(K1, K1)))

            K2q = np.cross(R, K1) + K2
            K2qt = np.cross(R, K1t)
            R13 = 1 / (np.linalg.norm(Rs - R) ** 3)
            R3 = R13 + 1 / (np.linalg.norm(R) ** 3)

            # K1N = K1/np.linalg.norm(K1)
            # P1.append(K1N)
            # if np.linalg.norm(K2) != 0:
            #     K2N = K2/np.linalg.norm(K2)
            # else:
            #     K2N = [0, 0, 0]
            # P2.append(K2N)

            # P1.append(np.linalg.norm(v))
            # P2.append(np.linalg.norm(K1))
            # P3.append(np.linalg.norm(K2))

            # print(np.dot(K1N, K2N))
            # print(np.dot(np.cross(R, K1), v + b/(m*j-b**2)*K2))
            # print(np.linalg.norm(v + b/(m*j-b**2)*K2))
            # print(np.dot(v + b/(m*j-b**2)*K2, v + b/(m*j-b**2)*K2))
            # print(np.dot(np.cross(R, K1), v + b/(m*j-b**2)*K2))
            # print(np.dot(v, v)+2*b*j/(m*j-b**2)**2*np.dot(K1, K2))
            # print(np.dot(K2t, K1))
            # print(np.dot(K1T, v+np.cross(w, R)))
            # print(np.dot(np.cross(R, K1T), w) + np.dot(K1T, v))

            # print(np.cross(Rt + b /(m*j - b**2) * K2, K1))

            # F = Rt + b /(m*j - b**2) * K2
            # print(np.dot(F, K1)/(np.sqrt(np.dot(F, F)) * np.sqrt(np.dot(K1, K1))))
            # print(np.sqrt(np.dot(F, F)))

            # print(np.dot(K2q, [1,0,0]))
            # print(np.linalg.norm(K2q))
            # print(np.dot(K2q, K2qt))

            # print('-'*20)
            # P1.append(np.linalg.norm(K1))



            y_n = y_n1
            t_n = t_n1

        return np.array(t), np.array(y), np.array(E), [1, 2, 3], [4, 5, 6]
        # return np.array(t), np.array(y), np.array(P1), np.array(P2), np.array(P3)  # , np.array(P4)