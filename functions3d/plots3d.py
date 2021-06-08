import numpy as np
import matplotlib.pyplot.cm as cm


def plot(fig, ux, uy, uz, lx, ly, lz, delt, delx, dely, delz, v0, t, num):
    #Start Points
    X = np.arange(0.0, lx + 2 * delx, delx)
    Y = np.arange(0.0, ly + 2 * dely, dely)
    Z = np.arange(0.0, lz + 2 * delz, delz)
    X, Y, Z = np.meshgrid(X, Y, Z)
    X = X.ravel()
    Y = Y.ravel()
    Z = Z.ravel()

    #Vectors
    U = ux[:uy.shape[0], :uz.shape[1], :ux.shape[2]]
    V = uy[:uy.shape[0], :uz.shape[1], :ux.shape[2]]
    W = uz[:uy.shape[0], :uz.shape[1], :ux.shape[2]]
    abs_u = np.sqrt(U * U + V * V + W * W) + 1e-8
    U = U / abs_u
    V = V / abs_u
    W = W / abs_u
    U = U.ravel()
    V = V.ravel()
    W = W.ravel()
    m = np.max(abs_u)
    print(np.max(abs_u / m))

    #Plot
    ax = fig.add_subplot(111,
                         projection='3d',
                         title = 't = {} [s]'.format(t * 0.1),
                         xlim=(-0.1, lx * 1.1),
                         ylim=(-0.1, ly * 1.1),
                         zlim=(-0.1, lz * 1.1))
    ax.set(xlabel='x',ylabel='y',zlabel='z')
    '''for x in range(int(lx / delx)):
        for y in range(int(ly / dely)):
            for z in range(int(lz / delz)):
                ax.quiver(x * delx,
                          y * dely,
                          z * delz,
                          ux[x, y, z] / 5,
                          uy[x, y, z] / 5,
                          uz[x, y, z] / 5,
                          (ux[x, y, z] ** 2 + uy[x, y, z] ** 2 + uz[x, y, z] ** 2) ** 0.5,
                          cmap='jet',
                          normalize=True)'''
    ax.quiver(X, Y, Z, U, V, W, cmap=cm.jet(abs_u / max(m, v0)))
#    ax.colorbar()
#    ax.show()
#    fig.savefig('./imgss/{:0=10}.png'.format(num))

