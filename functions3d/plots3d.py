import numpy as np
import matplotlib.pyplot as plt
"""Plot Air Velocity Field in Arrow Style.
"""


def plot(ux, uy, uz, delx, dely, delz, xx, yy, zz, ax):
    """Plot Air Velocity Field Arrow.

    Args:
        ux (np.ndarray): Air velocity in x dimension.
        uy (np.ndarray): Air velocity in y dimension.
        uz (np.ndarray): Air velocity in z dimension.
        delx (float): Infinitensimal displacement in x dimension.
        dely (float): Infinitensimal displacement in y dimension.
        delz (float): Infinitensimal displacement in z dimension.
        xx (np.ndarray): Cordinate array.
        yy (np.ndarray): Cordinate array.
        zz (np.ndarray): Cordinate array.
        ax (matplotlib.pyplot.axes): Ax object prepared in main function.

    Returns:
        im(matplotlib.pyplot.axes): Ploted ax object.
    """
    #Vectors
    u = ux[:uy.shape[0], :uz.shape[1], :ux.shape[2]]
    v = uy[:uy.shape[0], :uz.shape[1], :ux.shape[2]]
    w = uz[:uy.shape[0], :uz.shape[1], :ux.shape[2]]
    abs_u = np.sqrt(u * u + v * v + w * w) + 1e-8

    xx = xx.ravel()
    yy = yy.ravel()
    zz = zz.ravel()
    u = u.ravel()
    v = v.ravel()
    w = w.ravel()
    abs_u = abs_u.ravel()


    xx = xx * delx
    yy = yy * dely
    zz = zz * delz

    c = np.concatenate((abs_u, np.repeat(abs_u, 2)))
    c = plt.cm.jet(c)

    #Plot
    im = ax.quiver(xx, yy, zz, u, v, w, length=0.1, colors=c, cmap='jet', normalize=True)

    return im
