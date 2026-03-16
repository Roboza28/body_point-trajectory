from collections import namedtuple
import numpy as np


params_tuple = namedtuple('params' , 'm, j, b, A, yC2')


def angle_between(vector_1, vector_2):
    dot_product = np.dot(vector_1, vector_2)
    angle_between_vectors = np.arccos(dot_product / (np.linalg.norm(vector_1) * np.linalg.norm(vector_2)))
    return np.rad2deg(angle_between_vectors)


