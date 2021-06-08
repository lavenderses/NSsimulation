import numpy as np


def divergence(ux_ast, uy_ast, uz_ast, delt, delx, dely, delz):
    bef_ux = ux_ast[:-1, :, :]
    aft_ux = ux_ast[1:, :, :]

    bef_uy = uy_ast[:, :-1, :]
    aft_uy = uy_ast[:, 1:, :]

    bef_uz = uz_ast[:, :, :-1]
    aft_uz = uz_ast[:, :, 1:]

    div = ((aft_ux - bef_ux) / (dely * delz) +
           (aft_uy - bef_uy) / (delz * delx) +
           (aft_uz - bef_uz) / (delx * dely)) / delt

    return div


def jacobi(div, delt, delx, dely, delz, P, eps, w):
    P_calc = np.copy(P)
    cnst = (delx * dely * delz) ** 2 * 0.5 / (delx ** 2 + dely ** 2 + delz ** 2)

    bef_px = P[:-2, 1:-1, 1:-1]
    aft_px = P[2:, 1:-1, 1:-1]
    bef_py = P[1:-1, :-2, 1:-1]
    aft_py = P[1:-1, 2:, 1:-1]
    bef_pz = P[1:-1, 1:-1, :-2]
    aft_pz = P[1:-1, 1:-1, 2:]

    P_calc[1:-1, 1:-1, 1:-1] = ((aft_px + bef_px) / delx ** 2 +
                                (aft_py + bef_py) / dely ** 2 +
                                (aft_pz + bef_pz) / delz ** 2 - div[1:-1, 1:-1, 1:-1]) * cnst

    error = np.max(np.abs(P_calc[1:-1, 1:-1, 1:-1] - P[1:-1, 1:-1, 1:-1]))
    error = max(error, eps)

    return error, P_calc


def poisson(ux_ast, uy_ast, uz_ast, delt, delx, dely, delz, P, eps, w, count_max):
    error = 1
    count = 0
    div = divergence(ux_ast, uy_ast, uz_ast, delt, delx, dely, delz)

    while error > eps or count > count_max:
        error, P = jacobi(div, delt, delx, dely,  delz, P, eps, w)
        count += 1

    return P


def fix_u(ux_ast, uy_ast, uz_ast, P, delt, delx, dely, delz, p_rho):
    ux = np.zeros_like(ux_ast)
    uy = np.zeros_like(uy_ast)
    uz = np.zeros_like(uz_ast)

    bef_px = P[:-1, :, :]
    aft_px = P[1:, :, :]
    bef_py = P[:, :-1, :]
    aft_py = P[:, 1:, :]
    bef_pz = P[:, :, :-1]
    aft_pz = P[:, :, 1:]

    ux[1:-1, :, :] = ux_ast[1:-1, :, :] - (aft_px - bef_px) * delt / (delx * p_rho)
    uy[:, 1:-1, :] = uy_ast[:, 1:-1, :] - (aft_py - bef_py) * delt / (dely * p_rho)
    uz[:, :, 1:-1] = uz_ast[:, :, 1:-1] - (aft_pz - bef_pz) * delt / (delz * p_rho)

    return ux, uy, uz
