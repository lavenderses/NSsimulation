import numpy as np
from mpl_toolkits.mplot3d import Axes3D


def plot(ux, uy, uz, lx, ly, lz, delt, delx, dely, delz, xx, yy, zz, v0, t, num, ax):
    #Vectors
    u = ux[:uy.shape[0], :uz.shape[1], :ux.shape[2]]
    v = uy[:uy.shape[0], :uz.shape[1], :ux.shape[2]]
    w = uz[:uy.shape[0], :uz.shape[1], :ux.shape[2]]
    abs_u = np.sqrt(u * u + v * v + w * w) + 1e-8
    u = u / abs_u
    v = v / abs_u
    w = w / abs_u

    xx = xx.ravel()
    yy = yy.ravel()
    zz = zz.ravel()
    u = u.ravel()
    v = v.ravel()
    w = w.ravel()
    
    xx = xx * delx
    yy = yy * dely
    zz = zz * delz

    #Plot
    if ax is None:
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(16, 12))
        ax = fig.add_subplot(111,
                             projection='3d',
                             xlim=(-0.1, lx * 1.1),
                             ylim=(-0.1, ly * 1.1))

        ax.set_title('t = {} [s]'.format(t * 0.1))
        ax.colorbar()
        im = ax.quiver(xx, yy, zz, u, v, w,  cmap='jet', length=0.1)
        fig.colorbar(im)
        fig.savefig('./imgs/{:0=10}.png'.format(num))

        plt.clf()
        plt.cla()
        plt.close()
        return None

    else:
        ax.set_title('t = {} [s]'.format(t * 0.1))
        im = ax.quiver(xx, yy, zz, u, v, w,  cmap='jet', length=0.1)
        return im

