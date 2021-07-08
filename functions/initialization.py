import numpy as np
"""Initialize Each Variables.
"""


def force_and_first(ux, uy, delt, hs, v0, theta, w=None, rad_range=None):
    """Initialize velocity and boundary conditions.

    Args:
        ux (np.ndarray): Air velocity in x dimension.
        uy (np.ndarray)): Air velocity in y dimension.
        delt (float): Infinitensimal time.
        hs (list): Fan's cordinate.
        v0 (float): Fan's initial velocity.
        theta (float): Fan's elevation angle.
        w (float, optional): Fan's angular velocity. Defaults to None.
        rad_range (list, optional): Fan's angle range. Defaults to None.

    Returns:
        ux (np.ndarray): Air velocity in x dimension.
        uy (np.ndarray)): Air velocity in y dimension.
        theta (float): Fan's elevation angle.
        w (float, optional): Fan's angular velocity. Defaults to None.
    """
    # Boundary Conditions (v = 0 on the Walls)
    ux[:, [0, -1]] = 0
    uy[[0, -1], :] = 0
    ux[[0, 1, -2, -1], :] = 0
    ux[:, [0, 1, -2, -1]] = 0

    if w is None:
        rad_range = [-4 * np.pi, 4 * np.pi]
    else:
        if (theta >= rad_range[1] and w >0) or (theta <= rad_range[0] and w < 0):
            w *= -1
        theta += w * delt

    # Fan
    for h in hs:
        ux[h[0], h[1]] = v0 * np.cos(theta)
        uy[h[0], h[1]] = v0 * np.sin(theta)

    return ux, uy, theta, w
