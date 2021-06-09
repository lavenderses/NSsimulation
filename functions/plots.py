import numpy as np


def plot(ux, uy, lx, ly, delt, delx, dely, v0, t, num, ax):
    #Start Points
    x = [np.full(ux.shape[1], j * dely) for j in range(uy.shape[0])]
    y = [np.arange(0, ux.shape[1]) * delx] * (uy.shape[0])
    x = np.array(x)
    y = np.array(y)

    #Vectors
    u = ux[:uy.shape[0], :ux.shape[1]]
    v = uy[:uy.shape[0], :ux.shape[1]]
    abs_u = np.sqrt(u * u + v * v)
    u = u / abs_u
    v = v / abs_u

    #Plot
    if ax is None:
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(16, 12))
        ax = fig.add_subplot(111,
                             projection='3d',
                             xlim=(-0.1, lx * 1.1),
                             ylim=(-0.1, ly * 1.1))
    
        ax.set_title('t = {} [s]'.format(t * 0.1))
        im = ax.quiver(x, y, u, v, abs_u, cmap='jet', length=.5)
        fig.colorbar(im)
        fig.savefig('./imgs/{:0=10}.png'.format(num))

        plt.clf()
        plt.cla()
        plt.close()
        return None

    else:
        im = ax.quiver(x, y, u, v, abs_u, cmap='jet', length=0.1)
        return im
