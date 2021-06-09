import numpy as np

def force_and_first(ux, uy, uz, delt, hs, v0, theta, phai):
    #Boundary Conditions (v = 0 on the Walls)
    #ux (Floor and Ceiling and Walls)
    ux[:, [0, -1], [0, -1]] = 0
    ux[[0, -1], :, [0, -1]] = 0
    ux[[0, -1], [0, -1], :] = 0

    #uy (Floor and Ceiling and Walls)
    uy[:, [0, -1], [0, -1]] = 0
    uy[[0, -1], :, [0, -1]] = 0
    uy[[0, -1], [0, -1], :] = 0

    #uy (Floor and Ceiling and Walls)
    uz[:, [0, -1], [0, -1]] = 0
    uz[[0, -1], :, [0, -1]] = 0
    uz[[0, -1], [0, -1], :] = 0

    #Fan
    for h in hs:
        ux[h[0], h[1], h[2]] = v0 * np.sin(theta) * np.cos(phai)
        uy[h[0], h[1], h[2]] = v0 * np.sin(theta) * np.sin(phai)
        uz[h[0], h[1], h[2]] = v0 * np.cos(theta)

    return ux, uy, uz