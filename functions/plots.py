import numpy as np
import matplotlib.cm as cm

def plot(ux, uy, Lx, Ly, delt, delx, dely, v0, plt, t):
    '''
    X = np.delete(ux, Lx + 2, 0)
    Y = np.delete(uy, Ly + 2, 1)
    X = np.flatten(X)
    Y = np.flatten(Y)'''
    for x in range(Lx + 2):
        for y in range(Ly + 2):
            v = np.array([ux[x, y], uy[x, y]])
            abs_u = np.linalg.norm(v)
            img = plt.quiver(x * delx,
                      y * dely,
                      ux[x, y] * delt,
                      uy[x, y] * delt,
                      color=cm.jet(abs_u / v0),
                      )
    plt.title('t = {} [s]'.format(t * 0.1))
    plt.colorbar(img)
    plt.show()

    return img
