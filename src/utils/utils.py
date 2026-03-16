from collections import namedtuple
import numpy as np


params_tuple = namedtuple('params' , 'm, j, b, A, yC2')


def angle_between(vector_1, vector_2):
    dot_product = np.dot(vector_1, vector_2)
    angle_between_vectors = np.arccos(dot_product / (np.linalg.norm(vector_1) * np.linalg.norm(vector_2)))
    return np.rad2deg(angle_between_vectors)


def derivative_vector(vector_act: np.ndarray, vector_old: np.ndarray, dt: float):
    vector_x = (vector_act[0] - vector_old[0]) / dt
    vector_y = (vector_act[1] - vector_old[1]) / dt
    vector_z = (vector_act[2] - vector_old[2]) / dt

    vector_new = np.array([vector_x, vector_y, vector_z])

    return vector_new
