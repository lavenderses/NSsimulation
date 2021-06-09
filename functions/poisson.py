import numpy as np

def divergence(ux_ast, uy_ast, delt, delx, dely):
    bef_ux = ux_ast[:-1, :]
    aft_ux = ux_ast[1:, :]

    bef_uy = uy_ast[:, :-1]
    aft_uy = uy_ast[:, 1:]

    div = ((aft_ux - bef_ux) / delx + (aft_uy - bef_uy) / dely) / delt

    return div


def jacobi(div, delt, delx, dely, p, eps):
    p_calc = np.copy(p)
    cnst = pow(delx * dely, 2) * 0.5 / (pow(delx, 2) + pow(dely, 2))

    bef_px = p[:-2, 1:-1]
    aft_px = p[2:, 1:-1]
    bef_py = p[1:-1, :-2]
    aft_py = p[1:-1, 2:]

    p_calc[1:-1, 1:-1] = ((aft_px + bef_px) / delx ** 2 +
                          (aft_py + bef_py) / dely ** 2 - div[1:-1, 1:-1]) * cnst
    error = np.max(np.abs(p_calc[1:-1, 1:-1] - p[1:-1, 1:-1]))
    error = max(error, eps)

    return error, p_calc


def poisson(ux_ast, uy_ast, delt, delx, dely, p, eps, w, count_max):
    error = 1
    count = 0
    div = divergence(ux_ast, uy_ast, delt, delx, dely)

    while error > eps or count > count_max:
        error, p = jacobi(div, delt, delx, dely, p, eps)
        count += 1

    return p


def fix_u(ux_ast, uy_ast, p, delt, delx, dely, p_rho):
    ux = np.zeros_like(ux_ast)
    uy = np.zeros_like(uy_ast)

    bef_px = p[:-1, :]
    aft_px = p[1:, :]
    bef_py = p[:, :-1]
    aft_py = p[:, 1:]

    ux[1:-1, :] = ux_ast[1:-1, :] - (aft_px - bef_px) * delt / (delx * p_rho)
    uy[:, 1:-1] = uy_ast[:, 1:-1] - (aft_py - bef_py) * delt / (dely * p_rho)

    return ux, uy
