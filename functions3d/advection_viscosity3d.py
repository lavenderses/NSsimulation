import numpy as np
"""Calculate Advection & Viscosity Member in Navier Stokes Equation.
"""


def advection_x(ux, uy, uz, delt, delx, dely, delz):
    """Advection item in x dimension.

    Args:
        ux (np.ndarray): Air velocity in x dimension.
        uy (np.ndarray)): Air velocity in y dimension.
        uz (np.ndarray)): Air velocity in z dimension.
        delt (float): Infinitensimal time.
        delx (float): Infinitensimal displacement in x dimension.
        dely (float): Infinitensimal displacement in y dimension.
        delz (float): Infinitensimal displacement in z dimension.

    Returns:
        ux_ast (np.ndarray): Air velocity in x dimension.
    """
    ux_ast = np.zeros_like(ux)
    for x in range(1, ux.shape[0] - 1):
        for y in range(1, ux.shape[1] - 1):
            for z in range(1, ux.shape[2] - 1):
                vx = ux[x, y, z]
                vy = (uy[x, y, z] + uy[x - 1, y, z] + uy[x, y + 1, z] + uy[x - 1, y + 1, z]) / 4
                vz = (uz[x, y, z] + uz[x - 1, y, z] + uz[x, y, z + 1] + uz[x - 1, y, z + 1]) / 4

                if vx >= 0:
                    if vy >= 0:
                        if vz >= 0:
                            ux_ast[x, y, z] = ux[x, y, z] - delt * (
                                        vx * (ux[x, y, z] - ux[x - 1, y, z]) / delx + 
                                        vy * (ux[x, y, z] - ux[x, y - 1, z]) / dely + 
                                        vz * (ux[x, y, z] - ux[x, y, z - 1]) / delz)
                        else:
                            ux_ast[x, y, z] = ux[x, y, z] - delt * (
                                        vx * (ux[x, y, z] - ux[x - 1, y, z]) / delx + 
                                        vy * (ux[x, y, z] - ux[x, y - 1, z]) / dely + 
                                        vz * (ux[x, y, z + 1] - ux[x, y, z]) / delz)

                    elif vy < 0:
                        if vz >= 0:
                            ux_ast[x, y, z] = ux[x, y, z] - delt * (
                                        vx * (ux[x, y, z] - ux[x - 1, y, z]) / delx + 
                                        vy * (ux[x, y + 1, z] - ux[x, y, z]) / dely + 
                                        vz * (ux[x, y ,z] - ux[x, y, z - 1]) / delz)
                        else:
                            ux_ast[x, y, z] = ux[x, y, z] - delt * (
                                        vx * (ux[x, y, z] - ux[x - 1, y, z]) / delx + 
                                        vy * (ux[x, y + 1, z] - ux[x, y, z]) / dely + 
                                        vz * (ux[x, y ,z + 1] - ux[x, y, z]) / delz)

                elif vx < 0:
                    if vy >= 0:
                        if vz >= 0:
                            ux_ast[x, y, z] = ux[x, y, z] - delt * (
                                        vx * (ux[x + 1, y, z] - ux[x, y, z]) / delx + 
                                        vy * (ux[x, y, z] - ux[x, y - 1, z]) / dely + 
                                        vz * (ux[x, y, z] - ux[x, y, z - 1]) / delz)
                        else:
                            ux_ast[x, y, z] = ux[x, y, z] - delt * (
                                        vx * (ux[x + 1, y, z] - ux[x, y, z]) / delx + 
                                        vy * (ux[x, y, z] - ux[x, y - 1, z]) / dely + 
                                        vz * (ux[x, y, z + 1] - ux[x, y, z]) / delz)

                    else:
                        if vz >= 0:
                            ux_ast[x, y, z] = ux[x, y, z] - delt * (
                                        vx * (ux[x + 1, y, z] - ux[x, y, z]) / delx + 
                                        vy * (ux[x, y + 1, z] - ux[x, y, z]) / dely + 
                                        vz * (ux[x, y, z] - ux[x, y, z - 1]) / delz)
                        else:
                            ux_ast[x, y, z] = ux[x, y, z] - delt * (
                                        vx * (ux[x + 1, y, z] - ux[x, y, z]) / delx + 
                                        vy * (ux[x, y + 1, z] - ux[x, y, z]) / dely + 
                                        vz * (ux[x, y, z + 1] - ux[x, y, z]) / delz)

    return ux_ast


def advection_y(ux, uy, uz, delt, delx, dely, delz):
    """Advection item in y dimension.

    Args:
        ux (np.ndarray): Air velocity in x dimension.
        uy (np.ndarray)): Air velocity in y dimension.
        uz (np.ndarray)): Air velocity in z dimension.
        delt (float): Infinitensimal time.
        delx (float): Infinitensimal displacement in x dimension.
        dely (float): Infinitensimal displacement in y dimension.
        delz (float): Infinitensimal displacement in z dimension.

    Returns:
        uy_ast (np.ndarray): Air velocity in y dimension.
    """
    uy_ast = np.zeros_like(uy)
    for x in range(1, uy.shape[0] - 1):
        for y in range(1, uy.shape[1] - 1):
            for z in range(1, uy.shape[2] - 1):
                vx = (ux[x, y, z] + ux[x, y - 1, z] + ux[x, y, z + 1] + ux[x, y - 1, z + 1]) / 4
                vy = uy[x, y, z]
                vz = (uz[x, y, z] + uz[x, y - 1, z] + uz[x + 1, y, z] + uz[x + 1, y - 1, z]) / 4

                if vx >= 0:
                    if vy >= 0:
                        if vz >= 0:
                            uy_ast[x, y, z] = uy[x, y, z] - delt * (
                                        vx * (uy[x, y, z] - uy[x - 1, y, z]) / delx + 
                                        vy * (uy[x, y, z] - uy[x, y - 1, z]) / dely + 
                                        vz * (uy[x, y, z] - uy[x, y, z - 1]) / delz)
                        else:
                            uy_ast[x, y, z] = uy[x, y, z] - delt * (
                                        vx * (uy[x, y, z] - uy[x - 1, y, z]) / delx + 
                                        vy * (uy[x, y, z] - uy[x, y - 1, z]) / dely + 
                                        vz * (uy[x, y, z + 1] - uy[x, y, z]) / delz)

                    elif vy < 0:
                        if vz >= 0:
                            uy_ast[x, y, z] = uy[x, y, z] - delt * (
                                        vx * (uy[x, y, z] - uy[x - 1, y, z]) / delx + 
                                        vy * (uy[x, y + 1, z] - uy[x, y, z]) / dely + 
                                        vz * (uy[x, y ,z] - uy[x, y, z - 1]) / delz)
                        else:
                            uy_ast[x, y, z] = uy[x, y, z] - delt * (
                                        vx * (uy[x, y, z] - uy[x - 1, y, z]) / delx + 
                                        vy * (uy[x, y + 1, z] - uy[x, y, z]) / dely + 
                                        vz * (uy[x, y ,z + 1] - uy[x, y, z]) / delz)

                elif vx < 0:
                    if vy >= 0:
                        if vz >= 0:
                            uy_ast[x, y, z] = uy[x, y, z] - delt * (
                                        vx * (uy[x + 1, y, z] - uy[x, y, z]) / delx + 
                                        vy * (uy[x, y, z] - uy[x, y - 1, z]) / dely + 
                                        vz * (uy[x, y, z] - uy[x, y, z - 1]) / delz)
                        else:
                            uy_ast[x, y, z] = uy[x, y, z] - delt * (
                                        vx * (uy[x + 1, y, z] - uy[x, y, z]) / delx + 
                                        vy * (uy[x, y, z] - uy[x, y - 1, z]) / dely + 
                                        vz * (uy[x, y, z + 1] - uy[x, y, z]) / delz)

                    else:
                        if vz >= 0:
                            uy_ast[x, y, z] = uy[x, y, z] - delt * (
                                        vx * (uy[x + 1, y, z] - uy[x, y, z]) / delx + 
                                        vy * (uy[x, y + 1, z] - uy[x, y, z]) / dely + 
                                        vz * (uy[x, y, z] - uy[x, y, z - 1]) / delz)
                        else:
                            uy_ast[x, y, z] = uy[x, y, z] - delt * (
                                        vx * (uy[x + 1, y, z] - uy[x, y, z]) / delx + 
                                        vy * (uy[x, y + 1, z] - uy[x, y, z]) / dely + 
                                        vz * (uy[x, y, z + 1] - uy[x, y, z]) / delz)

    return uy_ast



def advection_z(ux, uy, uz, delt, delx, dely, delz):
    """Advection item in x dimension.

    Args:
        ux (np.ndarray): Air velocity in x dimension.
        uy (np.ndarray)): Air velocity in y dimension.
        uz (np.ndarray)): Air velocity in z dimension.
        delt (float): Infinitensimal time.
        delx (float): Infinitensimal displacement in x dimension.
        dely (float): Infinitensimal displacement in y dimension.
        delz (float): Infinitensimal displacement in z dimension.

    Returns:
        uz_ast(np.ndarray): Air velocity in z dimension.
    """
    uz_ast = np.zeros_like(uz)
    for x in range(1, uz.shape[0] - 1):
        for y in range(1, uz.shape[1] - 1):
            for z in range(1, uz.shape[2] - 1):
                vx = (ux[x, y, z] + ux[x, y, z - 1] + ux[x, y + 1, z] + ux[x, y + 1, z - 1]) / 4
                vy = (uy[x, y, z] + uy[x, y, z - 1] + uy[x + 1, y, z] + uy[x + 1, y, z - 1]) / 4
                vz = uz[x, y, z]

                if vx >= 0:
                    if vy >= 0:
                        if vz >= 0:
                            uz_ast[x, y, z] = uz[x, y, z] - delt * (
                                        vx * (uz[x, y, z] - uz[x - 1, y, z]) / delx + 
                                        vy * (uz[x, y, z] - uz[x, y - 1, z]) / dely + 
                                        vz * (uz[x, y, z] - uz[x, y, z - 1]) / delz)
                        else:
                            uz_ast[x, y, z] = uz[x, y, z] - delt * (
                                        vx * (uz[x, y, z] - uz[x - 1, y, z]) / delx + 
                                        vy * (uz[x, y, z] - uz[x, y - 1, z]) / dely + 
                                        vz * (uz[x, y, z + 1] - uz[x, y, z]) / delz)

                    elif vy < 0:
                        if vz >= 0:
                            uz_ast[x, y, z] = uz[x, y, z] - delt * (
                                        vx * (uz[x, y, z] - uz[x - 1, y, z]) / delx + 
                                        vy * (uz[x, y + 1, z] - uz[x, y, z]) / dely + 
                                        vz * (uz[x, y ,z] - uz[x, y, z - 1]) / delz)
                        else:
                            uz_ast[x, y, z] = uz[x, y, z] - delt * (
                                        vx * (uz[x, y, z] - uz[x - 1, y, z]) / delx + 
                                        vy * (uz[x, y + 1, z] - uz[x, y, z]) / dely + 
                                        vz * (uz[x, y ,z + 1] - uz[x, y, z]) / delz)

                elif vx < 0:
                    if vy >= 0:
                        if vz >= 0:
                            uz_ast[x, y, z] = uz[x, y, z] - delt * (
                                        vx * (uz[x + 1, y, z] - uz[x, y, z]) / delx + 
                                        vy * (uz[x, y, z] - uz[x, y - 1, z]) / dely + 
                                        vz * (uz[x, y, z] - uz[x, y, z - 1]) / delz)
                        else:
                            uz_ast[x, y, z] = uz[x, y, z] - delt * (
                                        vx * (uz[x + 1, y, z] - uz[x, y, z]) / delx + 
                                        vy * (uz[x, y, z] - uz[x, y - 1, z]) / dely + 
                                        vz * (uz[x, y, z + 1] - uz[x, y, z]) / delz)

                    else:
                        if vz >= 0:
                            uz_ast[x, y, z] = uz[x, y, z] - delt * (
                                        vx * (uz[x + 1, y, z] - uz[x, y, z]) / delx + 
                                        vy * (uz[x, y + 1, z] - uz[x, y, z]) / dely + 
                                        vz * (uz[x, y, z] - uz[x, y, z - 1]) / delz)
                        else:
                            uz_ast[x, y, z] = uz[x, y, z] - delt * (
                                        vx * (uz[x + 1, y, z] - uz[x, y, z]) / delx + 
                                        vy * (uz[x, y + 1, z] - uz[x, y, z]) / dely + 
                                        vz * (uz[x, y, z + 1] - uz[x, y, z]) / delz)

    return uz_ast


def viscosity(ux, uy, uz, ux_ast, uy_ast, uz_ast, delx, dely, delz, mu, p_rho):
    """Viscosity item in x dimension.

    Args:
        ux (np.ndarray): Air velocity in x dimension.
        uy (np.ndarray)): Air velocity in y dimension.
        uz (np.ndarray)): Air velocity in z dimension.
        ux_ast (np.ndarray): Air velocity after advection in x dimension.
        uy_ast (np.ndarray): Air velocity after advection in y dimension.
        uz_ast (np.ndarray): Air velocity after advection in z dimension.
        delt (float): Infinitensimal time.
        delx (float): Infinitensimal displacement in x dimension.
        dely (float): Infinitensimal displacement in y dimension.
        delz (float): Infinitensimal displacement in z dimension.
        MU (float): Viscosity value.
        p_rho (float): Density field.

    Returns:
        ux_ast (np.ndarray): Air velocity in x dimension.
        uy_ast (np.ndarray): Air velocity in y dimension.
        uz_ast (np.ndarray): Air velocity in z dimension.
    """
    nu = mu / p_rho
    #d2ux/dx2
    bef_ux_x = ux[:-2, 1:-1, 1:-1]
    mid_ux_x = ux[1:-1, 1:-1, 1:-1]
    aft_ux_x = ux[2:, 1:-1, 1:-1]

    #d2ux/dy2
    bef_ux_y = ux[1:-1, :-2, 1:-1]
    mid_ux_y = ux[1:-1, 1:-1, 1:-1]
    aft_ux_y = ux[1:-1, 2:, 1:-1]

    #d2ux/dz2
    bef_ux_z = ux[1:-1, 1:-1, :-2]
    mid_ux_z = ux[1:-1, 1:-1, 1:-1]
    aft_ux_z = ux[1:-1, 1:-1, 2:]

    #d2uy/dx2
    bef_uy_x = uy[:-2, 1:-1, 1:-1]
    mid_uy_x = uy[1:-1, 1:-1, 1:-1]
    aft_uy_x = uy[2:, 1:-1, 1:-1]

    #d2uy/dy2
    bef_uy_y = uy[1:-1, :-2, 1:-1]
    mid_uy_y = uy[1:-1, 1:-1, 1:-1]
    aft_uy_y = uy[1:-1, 2:, 1:-1]

    #d2ux/dz2
    bef_uy_z = uy[1:-1, 1:-1, :-2]
    mid_uy_z = uy[1:-1, 1:-1, 1:-1]
    aft_uy_z = uy[1:-1, 1:-1, 2:]

    #d2uz/dx2
    bef_uz_x = uz[:-2, 1:-1, 1:-1]
    mid_uz_x = uz[1:-1, 1:-1, 1:-1]
    aft_uz_x = uz[2:, 1:-1, 1:-1]

    #d2uz/dy2
    bef_uz_y = uz[1:-1, :-2, 1:-1]
    mid_uz_y = uz[1:-1, 1:-1, 1:-1]
    aft_uz_y = uz[1:-1, 2:, 1:-1]
    
    #d2uz/dz2
    bef_uz_z = uz[1:-1, 1:-1, :-2]
    mid_uz_z = uz[1:-1, 1:-1, 1:-1]
    aft_uz_z = uz[1:-1, 1:-1, 2:]

    ux_ast[1:-1, 1:-1, 1:-1] += ((aft_ux_x - 2 * mid_ux_x + bef_ux_x) / delx ** 2 +
                                 (aft_ux_y - 2 * mid_ux_y + bef_ux_y) / dely ** 2 +
                                 (aft_ux_z - 2 * mid_ux_z + bef_ux_z) / delz ** 2) * nu

    uy_ast[1:-1, 1:-1, 1:-1] += ((aft_uy_x - 2 * mid_uy_x + bef_uy_x) / delx ** 2 +
                                 (aft_uy_y - 2 * mid_uy_y + bef_uy_y) / dely ** 2 +
                                 (aft_uy_z - 2 * mid_uy_z + bef_uy_z) / delz ** 2) * nu

    uz_ast[1:-1, 1:-1, 1:-1] += ((aft_uz_x - 2 * mid_uz_x + bef_uz_x) / delx ** 2 +
                                 (aft_uz_y - 2 * mid_uz_y + bef_uz_y) / dely ** 2 +
                                 (aft_uz_z - 2 * mid_uz_z + bef_uz_z) / delz ** 2) * nu


    return ux_ast, uy_ast, uz_ast
