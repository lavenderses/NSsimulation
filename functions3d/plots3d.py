import numpy as np
import matplotlib
from mpl_toolkits.mplot3d import Axes3D


def plot(ux, uy, uz, lx, ly, lz, delt, delx, dely, delz, xx, yy, zz, v0, t, ax):
    #Vectors
    u = ux[:uy.shape[0], :uz.shape[1], :ux.shape[2]]
    v = uy[:uy.shape[0], :uz.shape[1], :ux.shape[2]]
    w = uz[:uy.shape[0], :uz.shape[1], :ux.shape[2]]

    xx = xx.ravel()
    yy = yy.ravel()
    zz = zz.ravel()
    u = u.ravel()
    v = v.ravel()
    w = w.ravel()

    xx = xx * delx
    yy = yy * dely
    zz = zz * delz

    colors = np.concatenate([u, v, w])
    print(np.max(u), np.max(v), np.max(w))
    cmap = matplotlib.cm.bwr

    #Plot
    ax.set_title('t = {} [s]'.format(t))
    im = ax.quiver(xx, yy, zz, u, v, w,  color=cmap(colors), length=0.1)

    return im
