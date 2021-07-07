import numpy as np


def plot_particle(ux, uy, uz, ps, lx, ly, lz, delt, delx, dely, delz, v0, ax):
    p_map = np.copy(ps)
    p_map[:, 0] /= delx
    p_map[:, 1] /= dely
    p_map[:, 2] /= delz
    p_map = p_map.astype(np.int)
    max_x = uy.shape[0]
    max_y = uz.shape[1]
    max_z = ux.shape[2]
    p_map[:, 0] = np.clip(p_map[:, 0], 0, max_x - 1)
    p_map[:, 1] = np.clip(p_map[:, 1], 0, max_y - 1)
    p_map[:, 2] = np.clip(p_map[:, 2], 0, max_z - 1)

    for i in range(len(ps)):
        p = p_map[i]
        ps[i, 0] += ux[p[0], p[1], p[2]] * delt * 1000
        ps[i, 1] += uy[p[0], p[1], p[2]] * delt * 1000
        ps[i, 2] += uz[p[0], p[1], p[2]] * delt * 1000

    img = ax.scatter(ps[:, 0], ps[:, 1], ps[:, 2], 'o', color='red',)

    return img
