import numpy as np
"""Plot Air Velocity Field in Arrow Style.
"""

def plot(ux, uy, delx, dely, xx, yy, ax):
    """Plot Air Velocity Field Arrow.

    Args:
        ux (np.ndarray): Air velocity in x dimension.
        uy (np.ndarray): Air velocity in y dimension.
        delt (float): Infinitensimal time.
        delx (float): Infinitensimal displacement in x dimension.
        dely (float): Infinitensimal displacement in y dimension.
        xx (np.ndarray): Cordinate array.
        yy (np.ndarray): Cordinate array.
        ax (matplotlib.pyplot.axes): Ax object prepared in main function.

    Returns:
        im (matplotlib.pyplot.axes): Ploted ax object.
    """
    u = ux[:uy.shape[0], :ux.shape[1]]
    v = uy[:uy.shape[0], :ux.shape[1]]
    abs_u = np.sqrt(u * u + v * v) + 1e-8
    u = u / abs_u
    v = v / abs_u

    xx = xx.ravel()
    yy = yy.ravel()
    u = u.ravel()
    v = v.ravel()
    abs_u = abs_u.ravel()

    xx = xx * delx
    yy = yy * dely

    im = ax.quiver(xx, yy, u, v, abs_u, cmap='jet', scale=50)

    return im
