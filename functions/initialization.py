import numpy as np

def force_and_first(ux, uy, delt, hs, v0, theta):
    #Boundary Conditions (v = 0 on the Walls)
    #ux (Floor and Ceiling and Walls)
    ux[:, [0, -1]] = 0
    ux[[0, -1], :] = 0

    #uy (Floor and Ceiling and Walls)
    uy[:, [0, -1]] = 0
    uy[[0, -1], :] = 0

    #Fan
    for h in hs:
        ux[h[0], h[1]] = v0 * np.cos(theta)
        uy[h[0], h[1]] = v0 * np.sin(theta)

    return ux, uy