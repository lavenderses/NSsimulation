import numpy as np

def divergence(ux, uy, div, delt, delx, dely):
    for x in range(div.shape[0]):
        for y in range(div.shape[1]):
            div[x][y] = ((ux[x + 2][y + 1] - ux[x + 1][y + 1]) / delx +
                         (uy[x + 1][y + 2] - uy[x + 1][y + 1]) / dely
                         ) / delt

    return div


def SOR(ux, uy, div, delt, delx, dely, P, eps, w):
    for x in range(1, P.shape[0] - 1):
        for y in range(1, P.shape[1] - 1):
            P_calced = ((P[x + 1][y] + P[x - 1][y]) / np.power(delx, 2) +
                        (P[x][y + 1] + P[x][y - 1]) / np.power(dely, 2) - 
                        div[x][y]
                        ) * np.power(delx * dely, 2) * 0.5 / (np.power(delx, 2) + np.power(dely, 2))
            error = max(abs(P[x][y] - P_calced), eps)
            P[x][y] = (1 - w) * P[x][y] + w * P_calced

    return error, P


def Poisson(ux, uy, div, delt, delx, dely, P, eps, w, count_max):
    error = 1
    count = 0
    div = divergence(ux, uy, div, delt, delx, dely)
    while error > eps or count > count_max:
        error, P = SOR(ux, uy, div, delt, delx, dely, P, eps, w)
        count += 1

    return P


def fix_u(ux, uy, P, delt, delx, dely, p_rho):
    for x in range(P.shape[0] - 1):
        for y in range(P.shape[1] - 1):
            ux[x + 2][y + 1] -= (P[x + 1][y] - P[x][y]) * delt / delx
            uy[x + 1][y + 2] -= (P[x][y + 1] - P[x][y]) * delt / dely

    return ux, uy