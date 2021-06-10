import numpy as np


def plot(ux, uy, lx, ly, delt, delx, dely, xx, yy, v0, t, num, ax, dirname=''):

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
    
    #Plot
    if ax is None:
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(16, 12))
        ax = fig.add_subplot(111,
                             xlim=(-0.1, lx * 1.1),
                             ylim=(-0.1, ly * 1.1))
        ax.set_title('t = {} [s]'.format(t * 0.1))
        im = ax.quiver(xx, yy, u, v, abs_u, cmap='jet')
        fig.colorbar(im)
        fig.savefig('{}/{:0=10}.png'.format(dirname, num))

        plt.clf()
        plt.cla()
        plt.close()
        return None

    else:
        im = ax.quiver(xx, yy, u, v, abs_u, cmap='jet')
        return im
