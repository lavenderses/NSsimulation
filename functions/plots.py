import numpy as np


def plot(ux, uy, lx, ly, delt, delx, dely, xx, yy, v0, num, ax):

    #Vectors
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

    #Plot
    ax.set_title('t = {} [s]'.format(num * delt))
    im = ax.quiver(xx, yy, u, v, abs_u, cmap='jet')

    return im
