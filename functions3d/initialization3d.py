import numpy as np
"""Initialize Each Variables.
"""


def force_and_first(ux, uy, uz, hs, v0, theta, phai):
    """Initialize velocity and boundary conditions.

    Args:
        ux (np.ndarray): Air velocity in x dimension.
        uy (np.ndarray)): Air velocity in y dimension.
        uz (np.ndarray)): Air velocity in z dimension.
        hs (list): Fan's cordinate.
        v0 (float): Fan's initial velocity.
        theta (float): Fan's elevation angle.
        phai (float): Fan's azimuth angle.

    Returns:
        ux (np.ndarray): Air velocity in x dimension.
        uy (np.ndarray)): Air velocity in y dimension.
        uz (np.ndarray)): Air velocity in z dimension.
    """
    # Boundary Conditions (v = 0 on the Walls)
    # ux (Floor and Ceiling and Walls)
    ux[[0, 1, -2, -1], :, :] = 0
    uy[:, [0, 1, -2, -1], :] = 0
    uz[:, :, [0, 1, -2, -1]] = 0

    ux[:, [0, -1], [0, -1]] = 0
    uy[[0, -1], [0, -1], :] = 0
    uz[[0, -1], [0, -1], :] = 0

    # Fan
    for h in hs:
        ux[h[0], h[1], h[2]] = v0 * np.sin(theta) * np.cos(phai)
        uy[h[0], h[1], h[2]] = v0 * np.sin(theta) * np.sin(phai)
        uz[h[0], h[1], h[2]] = v0 * np.cos(theta)

    return ux, uy, uz
