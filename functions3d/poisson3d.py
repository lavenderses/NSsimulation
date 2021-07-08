import numpy as np
"""Solve Poisson Equation.
"""


def divergence(ux_ast, uy_ast, uz_ast, delt, delx, dely, delz):
    """Calculate air velocity divergence.

    Args:
        ux_ast (np.ndarray): Air velocity after Advection & Viscosity in x dimension.
        uy_ast (np.ndarray): Air velocity after Advection & Viscosity in y dimension.
        uz_ast (np.ndarray): Air velocity after Advection & Viscosity in z dimension.
        delt (float): Infinitensimal time.
        delx (float): Infinitensimal displacement in x dimension.
        dely (float): Infinitensimal displacement in y dimension.
        delz (float): Infinitensimal displacement in z dimension.

    Returns:
        div (np.ndarray): Divergence.
    """
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


def jacobi(div, delx, dely, delz, P, eps):
    """Jacovi function judge whether pressure is converged.

    Args:
        div (np.ndarray): Divergence.
        delt (float): Infinitensimal time.
        delx (float): Infinitensimal displacement in x dimension.
        dely (float): Infinitensimal displacement in y dimension.
        delz (float): Infinitensimal displacement in z dimension.
        P (np.ndarray): Pressure.
        eps (float): Infinitensimal value for Jaccobi function.

    Returns:
        error (np.float): Max difference between Before & After Pressure.
        p_calc (np.ndarray): Culclated Pressure.
    """
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


def poisson(ux_ast, uy_ast, uz_ast, delt, delx, dely, delz, p, eps, count_max):
    """Solve poisson equation.

    Args:
        ux_ast (np.ndarray): Air velocity after Advection & Viscosity in x dimension.
        uy_ast (np.ndarray): Air velocity after Advection & Viscosity in y dimension.
        uz_ast (np.ndarray): Air velocity after Advection & Viscosity in z dimension.
        delt (float): Infinitensimal time.
        delx (float): Infinitensimal displacement in x dimension.
        dely (float): Infinitensimal displacement in y dimension.
        delz (float): Infinitensimal displacement in z dimension.
        p (np.ndarray): Pressure.
        eps (float): Infinitensimal value for Jaccobi function.
        count_max (int): Max iterations value.

    Returns:
        p (np.ndarray): Converged pessure.
    """
    error = 1
    count = 0
    div = divergence(ux_ast, uy_ast, uz_ast, delt, delx, dely, delz)

    while error > eps or count > count_max:
        error, p = jacobi(div, delx, dely,  delz, p, eps)
        count += 1

    return p


def fix_u(ux_ast, uy_ast, uz_ast, P, delt, delx, dely, delz, p_rho):
    """Fix air velocity by converged pressure.

    Args:
        ux_ast (np.ndarray): Air velocity after Advection & Viscosity in x dimension.
        uy_ast (np.ndarray): Air velocity after Advection & Viscosity in y dimension.
        uz_ast (np.ndarray): Air velocity after Advection & Viscosity in z dimension.
        p (np.ndarray): Converged pessure.
        delt (float): Infinitensimal time.
        delx (float): Infinitensimal displacement in x dimension.
        dely (float): Infinitensimal displacement in y dimension.
        delz (float): Infinitensimal displacement in z dimension.
        p_rho (float): Density field.

    Returns:
        ux (np.ndarray): The true (converged) air velocity in x dimension.
        uy (np.ndarray)): The true (converged) air velocity in y dimension.
        uz (np.ndarray)): The true (converged) air velocity in z dimension.
    """
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
