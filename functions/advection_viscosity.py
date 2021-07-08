import numpy as np
"""Calculate Advection & Viscosity Member in Navier Stokes Equation.
"""

def advection_x(ux, uy, delt, delx, dely):
    """Advection item in x dimension.

    Args:
        ux (np.ndarray): Air velocity in x dimension.
        uy (np.ndarray)): Air velocity in y dimension.
        delt (float): Infinitensimal time.
        delx (float): Infinitensimal displacement in x dimension.
        dely (float): Infinitensimal displacement in y dimension.

    Returns:
        ux (np.ndarray): Air velocity in x dimension.
        uy (np.ndarray)): Air velocity in y dimension.
    """
    ux_ast = np.zeros_like(ux)
    for x in range(1, ux.shape[0] - 1):
        for y in range(1, ux.shape[1] - 1):            
            vx = ux[x, y]
            vy = (uy[x, y] + uy[x - 1, y] + uy[x, y + 1] + uy[x - 1, y + 1]) / 4

            if vx >= 0 and vy >= 0:
                ux_ast[x, y] = ux[x, y] - delt * (
                        vx * (ux[x, y] - ux[x - 1, y]) / delx + 
                        vy * (ux[x, y] - ux[x, y - 1]) / dely
                        )

            elif vx >= 0 and vy < 0:
                ux_ast[x, y] = ux[x, y] - delt * (
                        vx * (ux[x, y] - ux[x - 1, y]) / delx + 
                        vy * (ux[x, y + 1] - ux[x, y]) / dely
                        )

            elif vx < 0 and vy >= 0:
                ux_ast[x, y] = ux[x, y] - delt * (
                        vx * (ux[x + 1, y] - ux[x, y]) / delx + 
                        vy * (ux[x, y] - ux[x, y - 1]) / dely
                        )

            else:
                ux_ast[x, y] = ux[x, y] - delt * (
                        vx * (ux[x + 1, y] - ux[x, y]) / delx + 
                        vy * (ux[x, y + 1] - ux[x, y]) / dely
                        )

    return ux_ast


def advection_y(ux, uy, delt, delx, dely):
    """Advection item in x dimension.

    Args:
        ux (np.ndarray): Air velocity in x dimension.
        uy (np.ndarray)): Air velocity in y dimension.
        delt (float): Infinitensimal time.
        delx (float): Infinitensimal displacement in x dimension.
        dely (float): Infinitensimal displacement in y dimension.

    Returns:
        ux (np.ndarray): Air velocity in x dimension.
        uy (np.ndarray)): Air velocity in y dimension.
    """
    uy_ast = np.zeros_like(uy)
    for x in range(1, uy.shape[0] - 1):
        for y in range(1, uy.shape[1] - 1):
            vx = (ux[x, y] + ux[x, y - 1] + ux[x + 1, y] + ux[x + 1, y - 1]) / 4
            vy = uy[x, y]

            if vx >= 0 and vy >= 0:
                uy_ast[x, y] = uy[x, y] - delt * (
                        vx * (uy[x, y] - uy[x - 1, y]) / delx + 
                        vy * (uy[x, y] - uy[x, y - 1]) / dely
                        )

            elif vx >= 0 and vy < 0:
                uy_ast[x, y] = uy[x, y] - delt * (
                        vx * (uy[x, y] - uy[x - 1, y]) / delx + 
                        vy * (uy[x, y + 1] - uy[x, y]) / dely
                        )

            elif vx < 0 and vy >= 0:
                uy_ast[x, y] = uy[x, y] - delt * (
                        vx * (uy[x + 1, y] - uy[x, y]) / delx + 
                        vy * (uy[x, y] - uy[x, y - 1]) / dely
                        )

            else:
                uy_ast[x, y] = uy[x, y] - delt * (
                        vx * (uy[x + 1, y] - uy[x, y]) / delx + 
                        vy * (uy[x, y + 1] - uy[x, y]) / dely
                        )

    return uy_ast


def viscosity(ux, uy, ux_ast, uy_ast, delx, dely, mu, p_rho):
    """Viscosity item in x dimension.

    Args:
        ux (np.ndarray): Air velocity in x dimension.
        uy (np.ndarray)): Air velocity in y dimension.
        ux_ast (np.ndarray): Air velocity after advection in x dimension.
        uy_ast (np.ndarray): Air velocity after advection in y dimension.
        delt (float): Infinitensimal time.
        delx (float): Infinitensimal displacement in x dimension.
        dely (float): Infinitensimal displacement in y dimension.
        MU (float): Viscosity value.
        p_rho (float): Density field.

    Returns:
        ux_ast (np.ndarray): Air velocity in x dimension.
        uy_ast (np.ndarray): Air velocity in y dimension.
    """
    nu = mu / p_rho
    #d2ux/dx2
    bef_ux_x = ux[:-2, 1:-1]
    mid_ux_x = ux[1:-1, 1:-1]
    aft_ux_x = ux[2:, 1:-1]

    #d2ux/dy2
    bef_ux_y = ux[1:-1, :-2]
    mid_ux_y = ux[1:-1, 1:-1]
    aft_ux_y = ux[1:-1, 2:]

    #d2uy/dx2
    bef_uy_x = uy[:-2, 1:-1]
    mid_uy_x = uy[1:-1, 1:-1]
    aft_uy_x = uy[2:, 1:-1]

    #d2uy/dy2
    bef_uy_y = uy[1:-1, :-2]
    mid_uy_y = uy[1:-1, 1:-1]
    aft_uy_y = uy[1:-1, 2:]

    ux_ast[1:-1, 1:-1] += ((aft_ux_x - 2 * mid_ux_x + bef_ux_x) / delx ** 2 +
                           (aft_ux_y - 2 * mid_ux_y + bef_ux_y) / dely ** 2) * nu

    uy_ast[1:-1, 1:-1] += ((aft_uy_x - 2 * mid_uy_x + bef_uy_x) / delx ** 2 +
                           (aft_uy_y - 2 * mid_uy_y + bef_uy_y) / dely ** 2) * nu


    return ux_ast, uy_ast
