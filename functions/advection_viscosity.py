import numpy as np


def advectionX(ux, uy, delt, delx, dely):
    ux_calced = np.zeros_like(ux)
    for x in range(1, ux.shape[0] - 1):
        for y in range(1, ux.shape[1] - 1):            
            vx = ux[x][y]
            vy = (uy[x][y] + uy[x - 1][y] + uy[x][y + 1] + uy[x - 1][y + 1]) / 4

            if vx >= 0 and vy >= 0:
                ux_calced[x][y] = ux[x][y] - delt * (
                        vx * (ux[x][y] - ux[x - 1][y]) / delx + 
                        vy * (ux[x][y] - ux[x][y - 1]) / dely
                    )

            elif vx >= 0 and vy < 0:
                ux_calced[x][y] = ux[x][y] - delt * (
                        vx * (ux[x][y] - ux[x - 1][y]) / delx + 
                        vy * (ux[x][y + 1] - ux[x][y]) / dely
                    )

            elif vx < 0 and vy >= 0:
                ux_calced[x][y] = ux[x][y] - delt * (
                        vx * (ux[x + 1][y] - ux[x][y]) / delx + 
                        vy * (ux[x][y] - ux[x][y - 1]) / dely
                    )

            else:
                ux_calced[x][y] = ux[x][y] - delt * (
                        vx * (ux[x + 1][y] - ux[x][y]) / delx + 
                        vy * (ux[x][y + 1] - ux[x][y]) / dely
                    )

    return ux_calced


def advectionY(ux, uy, delt, delx, dely):
    uy_calced = np.zeros_like(uy)
    for x in range(1, uy.shape[0] - 1):
        for y in range(1, uy.shape[1] - 1):
            vx = (ux[x][y] + ux[x][y - 1] + ux[x + 1][y] + ux[x + 1][y - 1]) / 4
            vy = uy[x][y]

            if vx >= 0 and vy >= 0:
                uy_calced[x][y] = uy[x][y] - delt * (
                        vx * (uy[x][y] - uy[x - 1][y]) / delx + 
                        vy * (uy[x][y] - uy[x][y - 1]) / dely
                    )

            elif vx >= 0 and vy < 0:
                uy_calced[x][y] = uy[x][y] - delt * (
                        vx * (uy[x][y] - uy[x - 1][y]) / delx + 
                        vy * (uy[x][y + 1] - uy[x][y]) / dely
                    )

            elif vx < 0 and vy >= 0:
                uy_calced[x][y] = uy[x][y] - delt * (
                        vx * (uy[x + 1][y] - uy[x][y]) / delx + 
                        vy * (uy[x][y] - uy[x][y - 1]) / dely
                    )

            else:
                uy_calced[x][y] = uy[x][y] - delt * (
                        vx * (uy[x + 1][y] - uy[x][y]) / delx + 
                        vy * (uy[x][y + 1] - uy[x][y]) / dely
                    )

    return uy_calced


def viscosity(ux, uy, mu, delt, delx, dely, p_rho):
    ux_calced = np.zeros_like(ux)
    uy_calced = np.zeros_like(uy)
    for x in range(1, ux.shape[0] - 1):
        for y in range(1, ux.shape[1] - 1):
            ux_calced[x][y] = ux[x][y] + (
                               (
                                ux[x + 1][y] - 2 * ux[x][y] + ux[x - 1][y]
                               ) / np.power(delx, 2) + 
                               (
                                ux[x][y + 1] - 2 * ux[x][y] + ux[x][y - 1]
                               ) / np.power(dely, 2)
                           ) * delt * mu / p_rho

    for x in range(1, uy.shape[0] - 1):
        for y in range(1, uy.shape[1] - 1):
            uy_calced[x][y] = uy[x][y] + (
                               (
                                uy[x + 1][y] - 2 * uy[x][y] + uy[x - 1][y]
                               ) / np.power(delx, 2) + 
                               (
                                uy[x][y + 1] - 2 * uy[x][y] + uy[x][y - 1]
                               ) / np.power(dely, 2)
                           ) * delt * mu / p_rho

    return ux_calced, uy_calced

