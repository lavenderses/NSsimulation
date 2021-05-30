import numpy as np
import matplotlib.pyplot as plt


def plot(ux, uy, Lx, Ly, delt, delx, dely, v0, t, num):
    #Start Points
    X = [np.full(ux.shape[1], j * dely) for j in range(uy.shape[0])]
    Y = [np.arange(0, ux.shape[1]) * delx] * (uy.shape[0])
    X = np.array(X)
    Y = np.array(Y)

    #Vectors
    U = ux[:uy.shape[0], :ux.shape[1]]
    V = uy[:uy.shape[0], :ux.shape[1]]
    abs_u = np.sqrt(U * U + V * V)
    U = U / abs_u
    V = V / abs_u

    #Plot
    plt.figure(figsize=(16, 12))
    plt.xlim(-0.1, Lx * delx * 1.1)
    plt.ylim(-0.1, Ly * delx * 1.1)

    plt.quiver(X, Y, U, V, abs_u, cmap='jet', minshaft=3, headwidth=12)
    plt.title('t = {} [s]'.format(t * 0.1))
    plt.colorbar()
    plt.savefig('./imgs/{:0=10}.png'.format(num))

    plt.clf()
    plt.cla()
    plt.close()
