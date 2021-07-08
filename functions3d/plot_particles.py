import numpy as np
"""Plot Evry Particle.
"""


def plot_particle(ux, uy, uz, lx, ly, lz, ps, delt, delx, dely, delz, ax):
    """Plot particles.

    Args:
        ux (np.ndarray): Air velocity in x dimension.
        uy (np.ndarray): Air velocity in y dimension.
        uz (np.ndarray): Air velocity in z dimension.
        lx (float): Room length in x dimenstion.
        ly (float): Room length in y dimenstion.
        lz (float): Room length in z dimenstion.
        ps (np.ndarray): Pressure.
        delt (float): Infinitensimal time.
        delx (float): Infinitensimal displacement in x dimension.
        dely (float): Infinitensimal displacement in y dimension.
        delz (float): Infinitensimal displacement in z dimension.
        ax (matplotlib.pyplot.axes): Ax object prepared in main function.

    Returns:
        img (matplotlib.pyplot.axes): Ploted ax object.
    """
    p_map = np.copy(ps)
    p_map[:, 0] /= delx
    p_map[:, 1] /= dely
    p_map[:, 2] /= delz
    p_map = p_map.astype(np.int)
    p_map[:, 0] = np.clip(p_map[:, 0], 0, uy.shape[0] - 1)
    p_map[:, 1] = np.clip(p_map[:, 1], 0, uz.shape[1] - 1)
    p_map[:, 2] = np.clip(p_map[:, 2], 0, ux.shape[2] - 1)

    for i in range(len(ps)):
        p = p_map[i]
        ps[i, 0] += ux[p[0], p[1], p[2]] * delt
        ps[i, 1] += uy[p[0], p[1], p[2]] * delt
        ps[i, 2] += uz[p[0], p[1], p[2]] * delt

        if (not  0 < ps[i, 0] < lx) or (not 0 < ps[i, 1] < ly) or (not 0 < ps[i, 2] < lz):
            ps[i, 0] = np.random.rand() * lx
            ps[i, 1] = np.random.rand() * ly
            ps[i, 2] = np.random.rand() * lz

    img = ax.scatter(ps[:, 0], ps[:, 1], ps[:, 2], 'o', color='red',)

    return img
