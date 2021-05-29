import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import Normalize

def plot(ux, uy, Lx, Ly, delt, delx, dely, v0, plt, t):
    X = ux[:Lx+2, :Ly+2]
    Y = ux[:Ly+2, :Ly+2]
    X = X.flatten()
    Y = Y.flatten()
    fig = plt.figure()
    #'''
    for x in range(Lx + 2):
        for y in range(Ly + 2):
            v = np.array([ux[x, y], uy[x, y]])
            abs_u = np.linalg.norm(v)
            img = plt.quiver(x * delx,
                      y * dely,
                      ux[x, y] * delt,
                      uy[x, y] * delt,
                      cmap='jet',
                      color=cm.jet(abs_u / v0),
                      norm=Normalize(0.0, v0)
                      )
    '''
    img = plt.quiver(X, Y, cmap='jet')
    '''
    plt.title('t = {} [s]'.format(t * 0.1))
    fig.colorbar(img)
    img.set_clim(0.0, v0)
    plt.plot()
    plt.show()

    return img
