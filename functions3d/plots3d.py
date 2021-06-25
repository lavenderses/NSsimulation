import numpy as np
import matplotlib.pyplot as plt


def plot(ux, uy, uz, lx, ly, lz, delt, delx, dely, delz, xx, yy, zz, v0, num, ax):
    #Vectors
    u = ux[:uy.shape[0], :uz.shape[1], :ux.shape[2]]
    v = uy[:uy.shape[0], :uz.shape[1], :ux.shape[2]]
    w = uz[:uy.shape[0], :uz.shape[1], :ux.shape[2]]
    abs_u = np.sqrt(u * u + v * v + w * w) + 1e-8
    '''
    u = u / abs_u
    v = v / abs_u
    w = w / abs_u'''

    xx = xx.ravel()
    yy = yy.ravel()
    zz = zz.ravel()
    u = u.ravel()
    v = v.ravel()
    w = w.ravel()
    abs_u = abs_u.ravel()

    '''
    xx = xx[::8]
    yy = yy[::8]
    zz = zz[::8]
    u = u[::8]
    v = v[::8]
    w = w[::8]'''

    xx = xx * delx
    yy = yy * dely
    zz = zz * delz

    c = np.concatenate((abs_u, np.repeat(abs_u, 2)))
    c = plt.cm.jet(c)

    #Plot
    ax.set_title('t = {:.2f} [s]'.format(num * delt))
    im = ax.quiver(xx, yy, zz, u, v, w, length=0.1, colors=c, cmap='jet', normalize=True)

    return im
