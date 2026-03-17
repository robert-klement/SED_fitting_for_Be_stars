import numpy as np


def calc_n(alpha_IR, lam):

    beta = 1.5

    g0 = np.array([0.0952, 0.1001, 0.1097, 0.1250, 0.1470, 0.1761])
    g1 = np.array([0.0215, 0.0421, 0.0639, 0.0858, 0.1071, 0.1269])
    g2 = np.array([0.0145, 0.0130, 0.0111, 0.0090, 0.0068, 0.0046])
    b0 = np.array([2.2125, 1.6304, 1.1316, 0.6927, 0.2964, -0.0690])
    b1 = np.array([-1.5290, -1.3884, -1.2866, -1.2128, -1.1585, -1.1185])
    b2 = np.array([0.0563, 0.0413, 0.0305, 0.0226, 0.0169, 0.0126])

    g = np.exp(g0 + g1 * np.log(lam) + g2 * (np.log(lam))**2)
    b = np.exp(b0 + b1 * np.log(lam) + b2 * (np.log(lam))**2)

    u = (1 / (g + b)) * ((g1 + 2 * g2 * np.log(lam)) * g + (b1 + 2 * b2 * np.log(lam)) * b)

    n = (beta * (alpha_IR - 4) + 2 * u - 4) / (2 * alpha_IR - 8)

    return n
