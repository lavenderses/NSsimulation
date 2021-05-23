import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def plot(ux, uy, Lx, Ly, delt, delx, dely, plots, v0):
    plt.xlim(-0.1, Lx + 0.1)
    plt.ylim(-0.1, Ly + 0.1)
    for x in range(uy.shape[0]):
        for y in range(ux.shape[1]):
            abs_u = (np.sqrt(np.power(ux[x][y], 2) + 
                     np.power(uy[x][y], 2) +
                     0.01)
                     )
            img = plt.quiver(x * delx,
                      y * dely,
                      ux[x][y] * delt,
                      uy[x][y] * delt,
                      color=cm.jet(abs_u / v0),
                      )
    plt.show()
    plots.append([img])

    return plots
