import numpy as np
import plotly.graph_objs as go
from plotly.offline import plot
import matplotlib.pyplot as plt
from abc import ABC, abstractmethod


class Solver(ABC):
    def __init__(self, system_of_task, inkrement, t_start, t_end, t_delta, y0, params_task):
        self.system_of_task = system_of_task
        self.inkrement = inkrement
        self.t_start = t_start
        self.t_end = t_end
        self.t_delta = t_delta
        self.params_task = params_task
        self.y0 = y0
        self.fig = go.Figure()

        self.t = None
        self.y = None
        self.E = None

    def solve(self):
        self.t = np.array([self.t_start])
        # self.y = np.array([self.y0[0:3], self.y0[3:6], self.y0[6:9]])
        self.y = []
        self.E = []
        self.y.append([self.y0[0:3], self.y0[3:6], self.y0[6:9]])

        t_n = self.t_start
        y_n = self.y0

        energy_0 = self.calc_prom_values(y_n)
        while t_n < self.t_end:
            y_n1 = y_n + self.inkrement(self.system_of_task, y_n, self.t_delta, self.params_task)
            t_n1 = t_n + self.t_delta

            # self.y = np.append(self.y, np.array([y_n1[0:3], y_n1[3:6], y_n1[6:9]]))
            self.y.append([y_n1[0:3], y_n1[3:6], y_n1[6:9]])
            self.t = np.append(self.t, t_n1)

            energy_n1 = self.calc_prom_values(y_n1)

            if not np.allclose(energy_0, energy_n1, rtol=0.01, atol=0.0):
                raise ValueError('Энергия меняется более, чем на один процент!\nРекомендуется взять шаг меньше')

            y_n = y_n1
            t_n = t_n1


        self.y = np.array(self.y)
        # print()

    def create_trajectory_for_plotly_graph(self, df, name):
        self.fig.add_trace(go.Scatter3d(x=df[0], y=df[1], z=df[2],
                                        mode='lines',
                                        name=name,
                                        opacity=0.5,
                                        line=dict(
                                            color=['black'],
                                            colorscale='Viridis')))

    def create_object_for_plotly_graph(self, df, name, size, color):
        self.fig.add_trace(go.Scatter3d(x=[df[0]], y=[df[1]], z=[df[2]],
                                        mode='markers',
                                        name=name,
                                        marker=dict(
                                            size=size,
                                            color=[color],
                                            colorscale='Viridis')))

    @abstractmethod
    def create_plotly_graph(self):
        pass

    def create_graph_energy(self):
        plt.title("Зависимость полной энергии от времени")
        plt.xlabel('t')
        plt.ylabel('E')
        plt.grid(True)
        ax = plt.gca()  # получаем текущий объект axes
        ax.yaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))
        # ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.1e'))
        # ax.yaxis.get_offset_text().set_visible(False)
        # ax.ticklabel_format(useOffset=False, style='plain')
        plt.plot(self.t, self.E, 'k')
        plt.show()

        # plt.title("Зависимость полной энергии от времени")
        # plt.xlabel('t')
        # plt.ylabel('E (×10⁻¹³)')  # Указываем масштаб в подписи
        # # Масштабируем данные для отображения
        # scale = 1e13
        # plt.plot(self.t, self.E * scale, 'k')
        # # plt.show()

    @abstractmethod
    def calc_prom_values(self, y_n):
        pass


class B2T1(Solver):
    def calc_prom_values(self, y_n):
        Rs = self.params_task.yC2
        R  = np.array([y_n[0], y_n[1], y_n[2]])
        K1 = np.array([y_n[3], y_n[4], y_n[5]])
        K2 = np.array([y_n[6], y_n[7], y_n[8]])

        K1t = -self.params_task.A / np.linalg.norm(R) ** 3 * R - self.params_task.A / np.linalg.norm(R - Rs) ** 3 * (
                    R - Rs)
        K2t = -self.params_task.b / (self.params_task.m * self.params_task.j - self.params_task.b ** 2) * np.cross(K1,
                                                                                                                   K2)
        Rt = 1 / (self.params_task.m * self.params_task.j - self.params_task.b ** 2) * (
                    self.params_task.j * K1 - self.params_task.b * K2)

        v = 1 / (self.params_task.m * self.params_task.j - self.params_task.b ** 2) * (
                    self.params_task.j * K1 - self.params_task.b * K2)
        w = 1 / (self.params_task.b ** 2 - self.params_task.m * self.params_task.j) * (
                    self.params_task.b * K1 - self.params_task.m * K2)

        K = 1 / 2 * self.params_task.m * np.dot(v, v) + self.params_task.b * np.dot(v,
                                                                                    w) + 1 / 2 * self.params_task.j * np.dot(
            w, w)
        # U = -self.params_task.A * (1 / (np.linalg.norm(R)))
        U = -self.params_task.A / np.linalg.norm(R) - self.params_task.A / np.linalg.norm(R - Rs)
        E_act = K + U
        self.E.append(E_act)

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

        return E_act

    def create_plotly_graph(self):
        self.create_trajectory_for_plotly_graph(df=[self.y[:, 0, 0], self.y[:, 0, 1], self.y[:, 0, 2]],
                                                name="траектория цели")

        self.create_object_for_plotly_graph(df=[self.y[-1, 0, 0], self.y[-1, 0, 1], self.y[-1, 0, 2]],
                                            name="актуальное положение цели", size=10, color='red')
        self.create_object_for_plotly_graph(df=[0, 0, 0],
                                            name="массивный центр № 1", size=20, color='green')
        self.create_object_for_plotly_graph(df=self.params_task.yC2,
                                            name="массивный центр № 2", size=20, color='blue')


        plot(self.fig)
        # fig.write_html('../../data/trajectory_rk_2b_1t.html')


class B1T1(Solver):
    def calc_prom_values(self, y_n):
        R  = np.array([y_n[0], y_n[1], y_n[2]])
        K1 = np.array([y_n[3], y_n[4], y_n[5]])
        K2 = np.array([y_n[6], y_n[7], y_n[8]])

        # K1t = - q / np.linalg.norm(R)**3 * R - q / np.linalg.norm(R - Rs)**3
        # K2t = -b / (m * j - b ** 2) * np.cross(K1, K2)

        # e_theory = np.array([np.cos(gamma_theory), np.sin(gamma_theory), 0])
        # e_practice = np.array([np.cos(gamma_practice), np.sin(gamma_practice), 0])

        v = 1 / (self.params_task.m * self.params_task.j - self.params_task.b ** 2) * (
                    self.params_task.j * K1 - self.params_task.b * K2)
        w = 1 / (self.params_task.b ** 2 - self.params_task.m * self.params_task.j) * (
                    self.params_task.b * K1 - self.params_task.m * K2)

        K = 1 / 2 * self.params_task.m * np.dot(v, v) + self.params_task.b * np.dot(v, w) + 1 / 2 * self.params_task.j * np.dot(w, w)
        U = -self.params_task.A * (1 / (np.linalg.norm(R)))
        self.E.append(K + U)

        K1N = K1 / np.linalg.norm(K1)
        # P1.append(K1N)
        K2N = K2 / np.linalg.norm(K2)
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
        return K + U

    def create_plotly_graph(self):
        self.create_trajectory_for_plotly_graph(df=[self.y[:, 0, 0], self.y[:, 0, 1], self.y[:, 0, 2]],
                                                name="траектория цели")

        self.create_object_for_plotly_graph(df=[self.y[-1, 0, 0], self.y[-1, 0, 1], self.y[-1, 0, 2]],
                                            name="актуальное положение цели", size=10, color='red')
        self.create_object_for_plotly_graph(df=self.params_task.yC[0],
                                            name="массивный центр № 1", size=20, color='green')
        self.create_object_for_plotly_graph(df=self.params_task.yC[1],
                                            name="массивный центр № 2", size=20, color='blue')


        plot(self.fig)
        # fig.write_html('../../data/trajectory_rk_2b_1t.html')


class B2T2(Solver):
    def solve(self):
        self.t = np.array([self.t_start])
        # self.y = np.array([self.y0[0:3], self.y0[3:6], self.y0[6:9]])
        self.y = []
        self.E = []
        self.y.append([self.y0[0:3], self.y0[3:6], self.y0[6:9]])

        t_n = self.t_start
        y_n = self.y0

        energy_0 = self.calc_prom_values(y_n)
        while t_n < self.t_end:
            y_n1 = y_n + self.inkrement(self.system_of_task, y_n, self.t_delta, self.params_task)
            t_n1 = t_n + self.t_delta

            # self.y = np.append(self.y, np.array([y_n1[0:3], y_n1[3:6], y_n1[6:9]]))
            self.y.append([y_n1[0:3], y_n1[3:6], y_n1[6:9]])
            self.t = np.append(self.t, t_n1)

            energy_n1 = self.calc_prom_values(y_n1)

            if not np.allclose(energy_0, energy_n1, rtol=0.01, atol=0.0):
                raise ValueError('Энергия меняется более, чем на один процент!\nРекомендуется взять шаг меньше')

            y_n = y_n1
            t_n = t_n1


        self.y = np.array(self.y)

    def calc_prom_values(self, y_n):
        R  = np.array([y_n[0], y_n[1], y_n[2]])
        K1 = np.array([y_n[3], y_n[4], y_n[5]])
        K2 = np.array([y_n[6], y_n[7], y_n[8]])

        # K1t = - q / np.linalg.norm(R)**3 * R - q / np.linalg.norm(R - Rs)**3
        # K2t = -b / (m * j - b ** 2) * np.cross(K1, K2)

        # e_theory = np.array([np.cos(gamma_theory), np.sin(gamma_theory), 0])
        # e_practice = np.array([np.cos(gamma_practice), np.sin(gamma_practice), 0])

        v = 1 / (self.params_task.m * self.params_task.j - self.params_task.b ** 2) * (
                    self.params_task.j * K1 - self.params_task.b * K2)
        w = 1 / (self.params_task.b ** 2 - self.params_task.m * self.params_task.j) * (
                    self.params_task.b * K1 - self.params_task.m * K2)

        return 0

    def create_plotly_graph(self):
        self.create_trajectory_for_plotly_graph(df=[self.y[:, 0, 0], self.y[:, 0, 1], self.y[:, 0, 2]],
                                                name="траектория цели № 1")

        self.create_object_for_plotly_graph(df=[self.y[-1, 0, 0], self.y[-1, 0, 1], self.y[-1, 0, 2]],
                                            name="актуальное положение цели", size=10, color='red')
        self.create_object_for_plotly_graph(df=[0, 0, 0],
                                            name="массивный центр № 1", size=20, color='green')


        plot(self.fig)