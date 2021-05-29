import numpy as np

def divergence(ux_ast, uy_ast, delt, delx, dely):
    bef_ux = ux_ast[:-1, :]
    aft_ux = ux_ast[1:, :]

    bef_uy = uy_ast[:, :-1]
    aft_uy = uy_ast[:, 1:]

    div = ((aft_ux - bef_ux) / delx + (aft_uy - bef_uy) / dely) / delt

    return div


def Jacobi(div, delt, delx, dely, P, eps, w):
    P_calc = np.copy(P)
    cnst = pow(delx * dely, 2) * 0.5 / (pow(delx, 2) + pow(dely, 2))

    bef_Px = P[:-2, 1:-1]
    aft_Px = P[2:, 1:-1]
    bef_Py = P[1:-1, :-2]
    aft_Py = P[1:-1, 2:]

    P_calc[1:-1, 1:-1] = ((aft_Px + bef_Px) / delx ** 2 +
                          (aft_Py + bef_Py) / dely ** 2 - div[1:-1, 1:-1]) * cnst
    error = np.max(np.abs(P_calc[1:-1, 1:-1] - P[1:-1, 1:-1]))
    error = max(error, eps)

    return error, P_calc


def Poisson(ux_ast, uy_ast, delt, delx, dely, P, eps, w, count_max):
    error = 1
    count = 0
    div = divergence(ux_ast, uy_ast, delt, delx, dely)

    while error > eps or count > count_max:
        error, P = Jacobi(div, delt, delx, dely, P, eps, w)
        count += 1

    return P


def fix_u(ux_ast, uy_ast, P, delt, delx, dely, p_rho):
    ux = np.zeros_like(ux_ast)
    uy = np.zeros_like(uy_ast)

    bef_Px = P[:-1, :]
    aft_Px = P[1:, :]
    bef_Py = P[:, :-1]
    aft_Py = P[:, 1:]

    ux[1:-1, :] = ux_ast[1:-1, :] - (aft_Px - bef_Px) * delt / (delx * p_rho)
    uy[:, 1:-1] = uy_ast[:, 1:-1] - (aft_Py - bef_Py) * delt / (dely * p_rho)

    return ux, uy
